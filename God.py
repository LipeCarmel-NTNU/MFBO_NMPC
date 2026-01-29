import time
import subprocess
from pathlib import Path
from typing import Tuple

import torch

from gpytorch.mlls import ExactMarginalLogLikelihood
from botorch.fit import fit_gpytorch_mll
from botorch.models import SingleTaskGP
from botorch.models.transforms import Normalize, Standardize
from botorch.optim import optimize_acqf
from botorch.sampling.normal import SobolQMCNormalSampler
from botorch.utils.multi_objective.box_decompositions.non_dominated import FastNondominatedPartitioning

from botorch.acquisition.multi_objective.logei import qLogNoisyExpectedHypervolumeImprovement
from botorch.acquisition.acquisition import AcquisitionFunction

# Usa tua interface existente
from matlab_interface import send_theta, read_results, RESULTS_FILE


# =========================
# CONFIG (ajuste aqui)
# =========================
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
DTYPE = torch.double

N_INIT = 10          # pontos iniciais
N_ITER = 40          # iterações BO
Q_BATCH = 1          # q candidates por iteração (1 = sequencial)

POLL_S = 1.0         # polling do results.csv
WAIT_MATLAB_S = 2.0  # espera caso results.csv ainda não exista

# Aquisição: log-space
ELL_Z = 0.25         # kernel Z (largura)
GAMMA_TIME = 1.0     # peso da penalização de tempo
EPS = 1e-6           # estabilidade numérica
WZ = 1.0             # peso do termo z no log (opcional)
WT = 1.0             # peso do termo tempo no log (opcional)

# Bounds (teu theta) Mudei!
# theta = [f, theta_p, theta_m, q1..q3, r_u1..r_u3, r_du1..r_du3]
LB = torch.tensor([0.10, 0.0, 0.0] + [-3.0]*9, dtype=DTYPE, device=DEVICE)
UB = torch.tensor([0.60, 15.0, 7.0] + [ 3.0]*9, dtype=DTYPE, device=DEVICE)
# Se quiser tentar iniciar MATLAB automaticamente (Windows):
START_MATLAB = False
MATLAB_CMD = "matlab"  # ou caminho completo


# =========================
# Utils
# =========================
def snap_theta(X: torch.Tensor) -> torch.Tensor:
    """Clamp + arredonda variáveis inteiras (theta_p, theta_m)."""
    X = X.clone()
    X[..., 0] = torch.clamp(X[..., 0], LB[0], UB[0])  # f
    X[..., 1] = torch.round(torch.clamp(X[..., 1], LB[1], UB[1]))  # theta_p int
    X[..., 2] = torch.round(torch.clamp(X[..., 2], LB[2], UB[2]))  # theta_m int
    X[..., 3:] = torch.clamp(X[..., 3:], LB[3:], UB[3:])
    return X


def wait_for_matlab_ready():
    """Espera o results/results.csv existir (MATLAB cria no início)."""
    while True:
        if Path(RESULTS_FILE).exists():
            return
        print(f"[wait] esperando MATLAB criar {RESULTS_FILE} ...")
        time.sleep(WAIT_MATLAB_S)


def last_nrows() -> int:
    ts, *_ = read_results()
    return len(ts)


def wait_for_new_row(prev_rows: int, target_theta: torch.Tensor, timeout_s: float = 1e9) -> Tuple[float, float, float]:
    """
    Espera o MATLAB adicionar uma nova linha.
    Retorna (SSE, SSdU, runtime_s) da linha que bate com o theta enviado.
    """
    t0 = time.time()
    while True:
        if time.time() - t0 > timeout_s:
            raise TimeoutError("Timeout esperando results.csv atualizar.")

        ts, sse, ssd_u, cost, runtime, theta_mat = read_results()
        if len(ts) <= prev_rows:
            time.sleep(POLL_S)
            continue

        # pode ter mais de 1 linha nova; procura a que bate com o theta
        for k in range(prev_rows, len(ts)):
            th = torch.tensor(theta_mat[k], dtype=DTYPE, device=DEVICE)
            if torch.allclose(th, target_theta.view(-1), atol=1e-6, rtol=0):
                return float(sse[k]), float(ssd_u[k]), float(runtime[k])

        time.sleep(POLL_S)


def evaluate_theta(theta: torch.Tensor) -> Tuple[float, float, float]:
    """Envia theta -> MATLAB e espera SSE/SSdU/runtime."""
    theta = snap_theta(theta).view(-1)
    prev = last_nrows()
    send_theta(theta.detach().cpu().tolist())
    sse, ssd, rt = wait_for_new_row(prev, theta)
    return sse, ssd, rt


def build_models(X: torch.Tensor, Y_obj: torch.Tensor, Y_time: torch.Tensor):
    """2 GPs: (obj) e (runtime)."""
    bounds = torch.stack([LB, UB]).to(device=DEVICE, dtype=DTYPE)
    input_tf = Normalize(d=X.shape[-1], bounds=bounds)
    obj_tf = Standardize(m=2)
    time_tf = Standardize(m=1)

    model_obj = SingleTaskGP(X, Y_obj, input_transform=input_tf, outcome_transform=obj_tf)
    model_time = SingleTaskGP(X, Y_time, input_transform=input_tf, outcome_transform=time_tf)

    mll_obj = ExactMarginalLogLikelihood(model_obj.likelihood, model_obj)
    mll_time = ExactMarginalLogLikelihood(model_time.likelihood, model_time)

    fit_gpytorch_mll(mll_obj)
    fit_gpytorch_mll(mll_time)

    return model_obj, model_time


class MFMOLogAcq(AcquisitionFunction):
    """
    log_alpha = qLogNEHVI + WZ*log(az+eps) - WT*gamma*log(E[t]+eps)

    Onde az(f)=exp(- (1-f)^2 / (2*ell_z^2)).
    """
    def __init__(self, qlognehvi, model_time, ell_z=0.25, gamma=1.0, eps=1e-6, wz=1.0, wt=1.0):
        super().__init__(model=qlognehvi.model)
        self.qlognehvi = qlognehvi
        self.model_time = model_time
        self.ell_z = ell_z
        self.gamma = gamma
        self.eps = eps
        self.wz = wz
        self.wt = wt

    def forward(self, X: torch.Tensor) -> torch.Tensor:
        # X: batch x q x d

        # (1) termo multiobjetivo: já é LOG
        log_hv = self.qlognehvi(X)  # (batch)

        # (2) Kernel Z em log
        f = X[..., 0]  # (batch x q)
        az = torch.exp(-((1.0 - f) ** 2) / (2.0 * (self.ell_z ** 2)))  # (batch x q) em (0,1]
        az = az.mean(dim=-1)  # agrega em q -> (batch)
        log_az = torch.log(az + self.eps)

        # (3) penalização de tempo em log
        tmean = self.model_time.posterior(X).mean  # (batch x q x 1)
        tmean = tmean.mean(dim=-2).squeeze(-1)     # (batch)
        log_time_pen = -self.gamma * torch.log(tmean + self.eps)

        return log_hv + self.wz * log_az + self.wt * log_time_pen


def propose_candidate(X: torch.Tensor, Y_obj: torch.Tensor, Y_time: torch.Tensor, q: int = 1) -> torch.Tensor:
    model_obj, model_time = build_models(X, Y_obj, Y_time)

    # Queremos MINIMIZAR SSE e SSdU.
    # Então maximizamos Y = -[SSE, SSdU]
    Y = Y_obj.detach()

    # ref_point para maximização: deve ser "pior" que os pontos observados (mais baixo)
    y_min = Y.min(dim=0).values
    y_range = (Y.max(dim=0).values - y_min).clamp_min(1e-6)
    ref_point = (y_min - 0.2 * y_range).tolist()

    _ = FastNondominatedPartitioning(
        ref_point=torch.tensor(ref_point, dtype=DTYPE, device=DEVICE),
        Y=Y,
    )

    sampler = SobolQMCNormalSampler(sample_shape=torch.Size([128]))

    qlognehvi = qLogNoisyExpectedHypervolumeImprovement(
        model=model_obj,
        ref_point=ref_point,
        X_baseline=X,
        sampler=sampler,
        prune_baseline=True,
        cache_pending=True,
        max_iep=256,
    )

    acq = MFMOLogAcq(
        qlognehvi=qlognehvi,
        model_time=model_time,
        ell_z=ELL_Z,
        gamma=GAMMA_TIME,
        eps=EPS,
        wz=WZ,
        wt=WT,
    )

    bounds = torch.stack([LB, UB]).to(device=DEVICE, dtype=DTYPE)
    cand, _ = optimize_acqf(
        acq_function=acq,
        bounds=bounds,
        q=q,
        num_restarts=20,
        raw_samples=1024,
        options={"batch_limit": 5, "maxiter": 250},
    )

    return snap_theta(cand)


def maybe_start_matlab():
    if not START_MATLAB:
        return

    print("[info] tentando iniciar MATLAB automaticamente...")
    try:
        subprocess.Popen(
            [MATLAB_CMD, "-batch", "try, main_BO; catch ME, disp(getReport(ME)); end;"],
            cwd=str(Path(__file__).resolve().parent),
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            shell=False,
        )
        print("[info] MATLAB iniciado (modo batch).")
    except Exception as e:
        print(f"[warn] falhou iniciar MATLAB automaticamente: {e}")
        print("[warn] inicie o MATLAB manualmente e rode main_BO.")


def main():
    torch.set_default_dtype(DTYPE)
    maybe_start_matlab()
    wait_for_matlab_ready()

    # ===== Inicialização Sobol =====
    sobol = torch.quasirandom.SobolEngine(dimension=12, scramble=True)
    X = LB + (UB - LB) * sobol.draw(N_INIT).to(device=DEVICE, dtype=DTYPE)
    X = snap_theta(X)

    # Y_obj = -[SSE, SSdU] (para maximizar)
    Y_obj = torch.zeros(N_INIT, 2, dtype=DTYPE, device=DEVICE)
    Y_time = torch.zeros(N_INIT, 1, dtype=DTYPE, device=DEVICE)

    print(f"[init] avaliando {N_INIT} pontos iniciais...")
    for i in range(N_INIT):
        sse, ssd, rt = evaluate_theta(X[i])
        Y_obj[i, 0] = -sse
        Y_obj[i, 1] = -ssd
        Y_time[i, 0] = rt
        print(f"  init {i+1:02d}/{N_INIT}: SSE={sse:.4g}, SSdU={ssd:.4g}, rt={rt:.2f}s, f={X[i,0].item():.3f}")

    # ===== Loop BO =====
    print(f"[bo] iniciando {N_ITER} iterações...")
    for it in range(N_ITER):
        cand = propose_candidate(X, Y_obj, Y_time, q=Q_BATCH).squeeze(0)  # (d,)
        sse, ssd, rt = evaluate_theta(cand)

        X = torch.cat([X, cand.view(1, -1)], dim=0)
        Y_obj = torch.cat([Y_obj, torch.tensor([[-sse, -ssd]], dtype=DTYPE, device=DEVICE)], dim=0)
        Y_time = torch.cat([Y_time, torch.tensor([[rt]], dtype=DTYPE, device=DEVICE)], dim=0)

        print(f"  iter {it+1:03d}/{N_ITER}: SSE={sse:.4g}, SSdU={ssd:.4g}, rt={rt:.2f}s, f={cand[0].item():.3f}")

    # ===== Resumo =====
    SSE_all = (-Y_obj[:, 0]).detach().cpu()
    SSdU_all = (-Y_obj[:, 1]).detach().cpu()

    best_sse_idx = torch.argmin(SSE_all)
    best_ssdu_idx = torch.argmin(SSdU_all)

    print("\n[done] resumo:")
    print(f"  Melhor SSE:  {SSE_all[best_sse_idx].item():.6g}  (linha {best_sse_idx.item()+1})")
    print(f"  Melhor SSdU: {SSdU_all[best_ssdu_idx].item():.6g}  (linha {best_ssdu_idx.item()+1})")
    print("  Resultados completos: results/results.csv (gerenciado pelo MATLAB).")


if __name__ == "__main__":
    main()
