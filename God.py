import csv
import time
import subprocess
from pathlib import Path
from typing import Tuple, Optional, List, Dict

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


# ============================================================
# CONFIG (ajuste aqui)
# ============================================================
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
DTYPE = torch.double

N_INIT = 20          # pontos iniciais (DOE)
N_ITER = 40000          # iterações BO (após DOE)
Q_BATCH = 1          # q candidates por iteração (1 = sequencial)

POLL_S = 1.0         # polling do results.csv
WAIT_MATLAB_S = 2.0  # espera caso results.csv ainda não exista

# Aquisição: log-space
ELL_Z = 0.25         # kernel Z (largura)
GAMMA_TIME = 1.0     # peso da penalização de tempo
EPS = 1e-6           # estabilidade numérica
WZ = 1.0             # peso do termo z no log
WT = 1.0             # peso do termo tempo no log

# Bounds (teu theta)
# theta = [f, theta_p, theta_m, q1..q3, r_u1..r_u3, r_du1..r_du3]
LB = torch.tensor([0.01, 0.0, 0.0] + [-3.0] * 9, dtype=DTYPE, device=DEVICE)
UB = torch.tensor([1.00, 15.0, 7.0] + [3.0] * 9, dtype=DTYPE, device=DEVICE)
FIXED_THETA: Dict[int, float] = {3: 0.0, 6: -1000.0, 7: -1000.0, 8: -1000.0}
for _idx in (6, 7, 8):
    LB[_idx] = -1000.0
THETA_D = int(LB.numel())
OPT_IDXS = [i for i in range(THETA_D) if i not in FIXED_THETA]
LB_OPT = LB[OPT_IDXS]
UB_OPT = UB[OPT_IDXS]

# Se quiser tentar iniciar MATLAB automaticamente (Windows):
START_MATLAB = False
MATLAB_CMD = "matlab"  # ou caminho completo

# Log extra (além do results.csv do MATLAB)
BASE_DIR = Path(__file__).resolve().parent
TRACE_FILE = BASE_DIR / "results" / "opt_trace.csv"   # histórico “do Python” (DOE+OPT)
TRACE_FLUSH_EVERY = 1                                # escreve a cada nova avaliação


# ============================================================
# Utils
# ============================================================
def snap_theta(X: torch.Tensor) -> torch.Tensor:
    """Clamp + arredonda variáveis inteiras (theta_p, theta_m)."""
    X = X.clone()
    X[..., 0] = torch.clamp(X[..., 0], LB[0], UB[0])  # f
    X[..., 1] = torch.round(torch.clamp(X[..., 1], LB[1], UB[1]))  # theta_p int
    X[..., 2] = torch.round(torch.clamp(X[..., 2], LB[2], UB[2]))  # theta_m int
    X[..., 3:] = torch.clamp(X[..., 3:], LB[3:], UB[3:])
    for idx, val in FIXED_THETA.items():
        X[..., idx] = val
    return X


def opt_to_theta(X_opt: torch.Tensor) -> torch.Tensor:
    """Map optimization-space points to full 12D theta."""
    shp = X_opt.shape[:-1] + (THETA_D,)
    X_theta = torch.empty(shp, dtype=DTYPE, device=DEVICE)
    X_theta[..., OPT_IDXS] = X_opt
    for idx, val in FIXED_THETA.items():
        X_theta[..., idx] = val
    return snap_theta(X_theta)


def theta_to_opt(X_theta: torch.Tensor) -> torch.Tensor:
    """Project full 12D theta to optimization dimensions."""
    return snap_theta(X_theta)[..., OPT_IDXS]


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


def _theta_close(a: torch.Tensor, b: torch.Tensor, atol: float = 1e-9) -> bool:
    return torch.allclose(a.view(-1), b.view(-1), atol=atol, rtol=0.0)


def wait_for_new_row(prev_rows: int, target_theta: torch.Tensor, timeout_s: float = 1e9) -> Tuple[float, float, float, float]:
    """
    Espera o MATLAB adicionar uma nova linha.
    Retorna (SSE, SSdU, J, runtime_s) da linha que bate com o theta enviado.
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
            if _theta_close(th, target_theta):
                return float(sse[k]), float(ssd_u[k]), float(cost[k]), float(runtime[k])

        time.sleep(POLL_S)


def evaluate_theta(theta: torch.Tensor) -> Tuple[float, float, float, float]:
    """Envia theta -> MATLAB e espera SSE/SSdU/J/runtime."""
    theta = snap_theta(theta).view(-1)
    prev = last_nrows()
    send_theta(theta.detach().cpu().tolist())
    sse, ssd, jval, rt = wait_for_new_row(prev, theta)
    return sse, ssd, jval, rt


def load_history_from_results() -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, int]:
    """
    Lê o results.csv (via read_results) e reconstrói:
      X       : (n, d)
      Y_obj   : (n, 2)   = -[SSE, SSdU]
      Y_time  : (n, 1)   = runtime_s (positivo)
    Retorna também n (número de linhas).
    """
    ts, sse, ssd_u, cost, runtime, theta_mat = read_results()
    n = len(ts)
    if n == 0:
        d = THETA_D
        X = torch.empty((0, d), dtype=DTYPE, device=DEVICE)
        Y_obj = torch.empty((0, 2), dtype=DTYPE, device=DEVICE)
        Y_time = torch.empty((0, 1), dtype=DTYPE, device=DEVICE)
        return X, Y_obj, Y_time, 0

    X = torch.tensor(theta_mat, dtype=DTYPE, device=DEVICE)
    X = snap_theta(X)

    Y_obj = torch.stack(
        [-torch.tensor(sse, dtype=DTYPE, device=DEVICE),
         -torch.tensor(ssd_u, dtype=DTYPE, device=DEVICE)],
        dim=-1
    )
    rt = torch.tensor(runtime, dtype=DTYPE, device=DEVICE).clamp_min(EPS)
    Y_time = rt.view(-1, 1)
    return X, Y_obj, Y_time, n


def ensure_trace_header():
    TRACE_FILE.parent.mkdir(parents=True, exist_ok=True)
    if TRACE_FILE.exists():
        return
    with TRACE_FILE.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "row_idx", "phase", "doe_idx", "bo_idx",
            "timestamp_py", "acq_value",
            "ELL_Z", "GAMMA_TIME", "WZ", "WT",
            "t_floor",
            "SSE", "SSdU", "J", "runtime_s",
            *[f"theta_{i+1}" for i in range(THETA_D)],
        ])


def count_trace_rows() -> int:
    if not TRACE_FILE.exists():
        return 0
    with TRACE_FILE.open("r", newline="") as f:
        r = csv.reader(f)
        rows = list(r)
    return max(0, len(rows) - 1)  # exclui header


def append_trace_row(row: Dict):
    ensure_trace_header()
    with TRACE_FILE.open("a", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            row["row_idx"], row["phase"], row["doe_idx"], row["bo_idx"],
            row["timestamp_py"], row["acq_value"],
            row["ELL_Z"], row["GAMMA_TIME"], row["WZ"], row["WT"],
            row["t_floor"],
            row["SSE"], row["SSdU"], row["J"], row["runtime_s"],
            *row["theta_list"],
        ])


def sync_trace_with_existing_results(X: torch.Tensor, Y_obj: torch.Tensor, Y_time: torch.Tensor):
    """
    Se já existir results.csv com linhas antigas, garante que o opt_trace.csv tenha
    pelo menos essas mesmas linhas (para retomada limpa).
    Para linhas antigas, acq_value fica vazio.
    """
    ensure_trace_header()
    n_trace = count_trace_rows()
    n_res = X.shape[0]
    if n_trace >= n_res:
        return

    # Precisamos também de SSE/SSdU/J do CSV — read_results já tem.
    ts, sse, ssd_u, cost, runtime, theta_mat = read_results()

    for idx in range(n_trace, n_res):
        phase = "DOE" if idx < N_INIT else "OPT"
        doe_idx = idx + 1 if idx < N_INIT else ""
        bo_idx = (idx - N_INIT + 1) if idx >= N_INIT else ""
        row = {
            "row_idx": idx + 1,
            "phase": phase,
            "doe_idx": doe_idx,
            "bo_idx": bo_idx,
            "timestamp_py": "",          # não temos (foi gerado antes)
            "acq_value": "",             # não temos
            "ELL_Z": ELL_Z,
            "GAMMA_TIME": GAMMA_TIME,
            "WZ": WZ,
            "WT": WT,
            "t_floor": float(torch.tensor(runtime[:idx+1]).min().clamp_min(EPS).item()),
            "SSE": float(sse[idx]),
            "SSdU": float(ssd_u[idx]),
            "J": float(cost[idx]),
            "runtime_s": float(runtime[idx]),
            "theta_list": [float(v) for v in theta_mat[idx]],
        }
        append_trace_row(row)


# ============================================================
# Model + Acq (corrigidos para tempo)
# ============================================================
def build_models(X_opt: torch.Tensor, Y_obj: torch.Tensor, Y_time: torch.Tensor):
    """
    2 GPs:
      - model_obj: Y_obj = -[SSE, SSdU] (standardized)
      - model_time_log: log(runtime) (standardized)  -> garante positividade via lognormal
    """
    bounds = torch.stack([LB_OPT, UB_OPT]).to(device=DEVICE, dtype=DTYPE)
    input_tf = Normalize(d=X_opt.shape[-1], bounds=bounds)
    obj_tf = Standardize(m=2)
    time_tf = Standardize(m=1)

    # runtime positivo -> log
    Y_time_log = torch.log(Y_time.clamp_min(EPS))

    model_obj = SingleTaskGP(X_opt, Y_obj, input_transform=input_tf, outcome_transform=obj_tf)
    model_time_log = SingleTaskGP(X_opt, Y_time_log, input_transform=input_tf, outcome_transform=time_tf)

    mll_obj = ExactMarginalLogLikelihood(model_obj.likelihood, model_obj)
    mll_time = ExactMarginalLogLikelihood(model_time_log.likelihood, model_time_log)

    fit_gpytorch_mll(mll_obj)
    fit_gpytorch_mll(mll_time)

    return model_obj, model_time_log


class MFMOLogAcq(AcquisitionFunction):
    """
    log_alpha = qLogNEHVI + WZ*log(az+eps) - WT*gamma*log(E[t]+eps)

    Tempo:
      o GP é treinado em log(t). Assumindo lognormal:
        E[t] = exp(mu + 0.5*var)
    """
    def __init__(
        self,
        qlognehvi,
        model_time_log,
        ell_z=0.25,
        gamma=1.0,
        eps=1e-6,
        wz=1.0,
        wt=1.0,
        t_floor: float = 1e-6,
        max_log_et: float = 50.0,
    ):
        super().__init__(model=qlognehvi.model)
        self.qlognehvi = qlognehvi
        self.model_time_log = model_time_log
        self.ell_z = ell_z
        self.gamma = gamma
        self.eps = eps
        self.wz = wz
        self.wt = wt
        self.t_floor = float(t_floor)
        self.max_log_et = float(max_log_et)

    def forward(self, X: torch.Tensor) -> torch.Tensor:
        # (1) termo multiobjetivo: já é LOG
        log_hv = self.qlognehvi(X)  # (batch)
        log_hv = torch.nan_to_num(log_hv, neginf=-1e6, posinf=1e6)

        # (2) Kernel Z em log
        f = X[..., 0]  # (batch x q)
        az = torch.exp(-((1.0 - f) ** 2) / (2.0 * (self.ell_z ** 2)))  # (batch x q) em (0,1]
        az = az.mean(dim=-1)  # (batch)
        log_az = torch.log(az + self.eps)
        log_az = torch.nan_to_num(log_az, neginf=-1e6, posinf=1e6)

        # (3) penalização de tempo em log, com GP em log(t)
        post = self.model_time_log.posterior(X)
        mu = post.mean          # (batch x q x 1)
        var = post.variance     # (batch x q x 1)

        # log(E[t]) = mu + 0.5*var (lognormal)
        log_et = (mu + 0.5 * var).mean(dim=-2).squeeze(-1)  # (batch)
        log_et = torch.nan_to_num(log_et, neginf=-1e6, posinf=1e6).clamp_max(self.max_log_et)

        # E[t] (para aplicar piso observado)
        et = torch.exp(log_et).clamp_min(self.t_floor)

        log_time_pen = -self.gamma * torch.log(et + self.eps)
        log_time_pen = torch.nan_to_num(log_time_pen, neginf=-1e6, posinf=1e6)

        out = log_hv + self.wz * log_az + self.wt * log_time_pen
        return torch.nan_to_num(out, neginf=-1e6, posinf=1e6)


def propose_candidate(
    X: torch.Tensor,
    Y_obj: torch.Tensor,
    Y_time: torch.Tensor,
    q: int = 1
) -> Tuple[torch.Tensor, float, float, List[float]]:
    """
    Retorna:
      cand_opt: (q, d_opt) snapped
      acq_value: float (já avaliado no snapped)
      t_floor: float usado
      ref_point: list[float] usado no NEHVI
    """
    X_opt = theta_to_opt(X)
    model_obj, model_time_log = build_models(X_opt, Y_obj, Y_time)

    # ref_point para maximização: deve ser "pior" que os pontos observados (mais baixo)
    y_min = Y_obj.min(dim=0).values
    y_range = (Y_obj.max(dim=0).values - y_min).clamp_min(1e-6)
    ref_point = (y_min - 0.2 * y_range).tolist()

    _ = FastNondominatedPartitioning(
        ref_point=torch.tensor(ref_point, dtype=DTYPE, device=DEVICE),
        Y=Y_obj,
    )

    sampler = SobolQMCNormalSampler(sample_shape=torch.Size([128]))

    qlognehvi = qLogNoisyExpectedHypervolumeImprovement(
        model=model_obj,
        ref_point=ref_point,
        X_baseline=X_opt,
        sampler=sampler,
        prune_baseline=True,
        cache_pending=True,
        max_iep=256,
    )

    # piso = melhor tempo observado até agora
    t_floor = float(Y_time.min().clamp_min(EPS).item())

    acq = MFMOLogAcq(
        qlognehvi=qlognehvi,
        model_time_log=model_time_log,
        ell_z=ELL_Z,
        gamma=GAMMA_TIME,
        eps=EPS,
        wz=WZ,
        wt=WT,
        t_floor=t_floor,
    )

    bounds = torch.stack([LB_OPT, UB_OPT]).to(device=DEVICE, dtype=DTYPE)
    cand_opt, _ = optimize_acqf(
        acq_function=acq,
        bounds=bounds,
        q=q,
        num_restarts=20,
        raw_samples=1024,
        options={"batch_limit": 5, "maxiter": 250},
    )

    cand_opt = theta_to_opt(opt_to_theta(cand_opt))

    # avalia acq no snapped (pra log)
    with torch.no_grad():
        av = acq(cand_opt).detach().cpu().view(-1)[0].item()

    return cand_opt, float(av), t_floor, ref_point


# ============================================================
# DOE helper (para completar DOE se CSV tiver < N_INIT)
# ============================================================
def sobol_unique_points(n_new: int, existing_X: torch.Tensor, seed: int = 1234) -> torch.Tensor:
    """
    Gera pontos Sobol e evita duplicar os já existentes.
    """
    d = len(OPT_IDXS)
    sobol = torch.quasirandom.SobolEngine(dimension=d, scramble=True, seed=seed)

    pts = []
    tries = 0
    max_tries = 20000
    while len(pts) < n_new and tries < max_tries:
        tries += 1
        x_opt = LB_OPT + (UB_OPT - LB_OPT) * sobol.draw(1).to(device=DEVICE, dtype=DTYPE)
        x = theta_to_opt(opt_to_theta(x_opt)).view(1, -1)

        duplicate = False
        if existing_X.numel() > 0:
            # checa contra existing e contra os já selecionados
            for j in range(existing_X.shape[0]):
                if _theta_close(existing_X[j], x[0], atol=1e-9):
                    duplicate = True
                    break
        if not duplicate:
            for y in pts:
                if _theta_close(y[0], x[0], atol=1e-9):
                    duplicate = True
                    break

        if not duplicate:
            pts.append(x)

    if len(pts) < n_new:
        raise RuntimeError(f"Não consegui gerar {n_new} pontos Sobol únicos (gerou {len(pts)}).")

    return torch.cat(pts, dim=0)


def maybe_start_matlab():
    if not START_MATLAB:
        return
    print("[info] tentando iniciar MATLAB automaticamente...")
    try:
        subprocess.Popen(
            [MATLAB_CMD, "-batch", "try, main_BO; catch ME, disp(getReport(ME)); end;"],
            cwd=str(BASE_DIR),
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            shell=False,
        )
        print("[info] MATLAB iniciado (modo batch).")
    except Exception as e:
        print(f"[warn] falhou iniciar MATLAB automaticamente: {e}")
        print("[warn] inicie o MATLAB manualmente e rode main_BO.")


# ============================================================
# Main
# ============================================================
def main():
    torch.set_default_dtype(DTYPE)
    maybe_start_matlab()
    wait_for_matlab_ready()

    # 1) carrega histórico (DOE + OPT) do results.csv
    X, Y_obj, Y_time, n_rows = load_history_from_results()

    # 2) garante trace alinhado com o que já existe no results.csv
    if n_rows > 0:
        sync_trace_with_existing_results(X, Y_obj, Y_time)

    n_doe_done = min(n_rows, N_INIT)
    n_opt_done = max(0, n_rows - N_INIT)

    print(f"[resume] encontrado {n_rows} linhas em {Path(RESULTS_FILE)} -> DOE: {n_doe_done}/{N_INIT}, OPT: {n_opt_done}/{N_ITER}")

    # 3) completa DOE se necessário
    if n_rows < N_INIT:
        n_missing = N_INIT - n_rows
        print(f"[doe] faltam {n_missing} pontos para completar DOE (N_INIT={N_INIT}) ...")

        X_new_opt = sobol_unique_points(n_missing, existing_X=theta_to_opt(X), seed=1234)
        for i in range(n_missing):
            theta_opt = X_new_opt[i]
            theta = opt_to_theta(theta_opt)
            sse, ssd, jval, rt = evaluate_theta(theta)
            X = torch.cat([X, theta.view(1, -1)], dim=0)
            Y_obj = torch.cat([Y_obj, torch.tensor([[-sse, -ssd]], dtype=DTYPE, device=DEVICE)], dim=0)
            Y_time = torch.cat([Y_time, torch.tensor([[rt]], dtype=DTYPE, device=DEVICE)], dim=0)

            global_row = X.shape[0]
            row = {
                "row_idx": global_row,
                "phase": "DOE",
                "doe_idx": global_row,
                "bo_idx": "",
                "timestamp_py": int(time.time()),
                "acq_value": "",
                "ELL_Z": ELL_Z,
                "GAMMA_TIME": GAMMA_TIME,
                "WZ": WZ,
                "WT": WT,
                "t_floor": float(Y_time.min().clamp_min(EPS).item()),
                "SSE": float(sse),
                "SSdU": float(ssd),
                "J": float(jval),
                "runtime_s": float(rt),
                "theta_list": [float(v) for v in theta.view(-1).detach().cpu().tolist()],
            }
            append_trace_row(row)

            print(f"  DOE {global_row:02d}/{N_INIT}: SSE={sse:.4g}, SSdU={ssd:.4g}, rt={rt:.2f}s, f={theta[0].item():.3f}")

        n_opt_done = max(0, X.shape[0] - N_INIT)

    # 4) continua BO do ponto onde parou
    remaining = max(0, N_ITER - n_opt_done)
    if remaining == 0:
        print("[bo] nada a fazer: já completou todas as iterações.")
    else:
        print(f"[bo] continuando: faltam {remaining} iterações (de {N_ITER}) ...")

    for k in range(remaining):
        bo_idx = n_opt_done + k + 1  # 1-based dentro do BO
        cand_opt, acq_val, t_floor, ref_point = propose_candidate(X, Y_obj, Y_time, q=Q_BATCH)
        cand_opt = cand_opt.squeeze(0)  # (d_opt,)
        cand = opt_to_theta(cand_opt.view(1, -1)).squeeze(0)  # (d_full,)

        sse, ssd, jval, rt = evaluate_theta(cand)

        X = torch.cat([X, cand.view(1, -1)], dim=0)
        Y_obj = torch.cat([Y_obj, torch.tensor([[-sse, -ssd]], dtype=DTYPE, device=DEVICE)], dim=0)
        Y_time = torch.cat([Y_time, torch.tensor([[rt]], dtype=DTYPE, device=DEVICE)], dim=0)

        global_row = X.shape[0]
        row = {
            "row_idx": global_row,
            "phase": "OPT",
            "doe_idx": "",
            "bo_idx": bo_idx,
            "timestamp_py": int(time.time()),
            "acq_value": float(acq_val),
            "ELL_Z": ELL_Z,
            "GAMMA_TIME": GAMMA_TIME,
            "WZ": WZ,
            "WT": WT,
            "t_floor": float(t_floor),
            "SSE": float(sse),
            "SSdU": float(ssd),
            "J": float(jval),
            "runtime_s": float(rt),
            "theta_list": [float(v) for v in cand.view(-1).detach().cpu().tolist()],
        }
        append_trace_row(row)

        print(f"  OPT {bo_idx:03d}/{N_ITER}: SSE={sse:.4g}, SSdU={ssd:.4g}, rt={rt:.2f}s, f={cand[0].item():.3f}, acq={acq_val:.3g}")

    # 5) Resumo
    if X.shape[0] > 0:
        SSE_all = (-Y_obj[:, 0]).detach().cpu()
        SSdU_all = (-Y_obj[:, 1]).detach().cpu()

        best_sse_idx = torch.argmin(SSE_all)
        best_ssdu_idx = torch.argmin(SSdU_all)

        print("\n[done] resumo:")
        print(f"  Melhor SSE:  {SSE_all[best_sse_idx].item():.6g}  (linha {best_sse_idx.item()+1})")
        print(f"  Melhor SSdU: {SSdU_all[best_ssdu_idx].item():.6g}  (linha {best_ssdu_idx.item()+1})")
        print(f"  results.csv : {Path(RESULTS_FILE)}")
        print(f"  trace.csv   : {TRACE_FILE}")
    else:
        print("[done] nenhum dado em results.csv.")


if __name__ == "__main__":
    main()
