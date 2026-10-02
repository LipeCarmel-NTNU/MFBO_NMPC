"""A real-time abort becomes an imputed dominated row, in both phases.

    python tests/python/test_realtime.py

MATLAB is replaced by fakes of send_request, wait_for_result and
wait_for_matlab_ready. The baseline case avoids the phi fit, so no .mat files
are needed. The results go to a temporary folder.
"""
import os
import shutil
import sys
import tempfile
import time
from pathlib import Path

from _check import REPO, check, finish

WORK = Path(tempfile.mkdtemp(prefix="rt_py_"))
os.environ["MFBO_RESULTS_DIR"] = str(WORK)
sys.path.insert(0, str(REPO))
sys.argv = [sys.argv[0], "--case", "baseline"]

import run_config  # noqa: E402
from pipeline import driver  # noqa: E402
from pipeline import matlab_interface as mi  # noqa: E402


class Stop(Exception):
    pass


RT_ID = "NMPC:realtimeInfeasible"
MATLAB_SECONDS = 0.3


def measured_row(eval_id, phase, theta, sse):
    return {"eval_id": eval_id, "timestamp": f"20261002_0000{eval_id:02d}",
            "phase": phase, "phi_vintage": 0.0, "z": theta[0],
            "SSE_measured": sse, "SSdU_measured": 1e-3, "phi_SSE": 1.0, "phi_SSdU": 1.0,
            "SSE": sse, "SSdU": 1e-3, "J": sse + 10.0, "runtime_s": 100.0 + eval_id,
            "n_flag_not_one": 0.0, "phi_floored": 0.0, "wall_total_s": 120.0,
            "wall_cases_s": 110.0, "wall_phi_s": 0.01, "wall_build_s": 0.5,
            "wall_save_s": 0.2, "theta": theta}


def make_fakes(phase_key, phase_label, infeasible_ids, stop_after=None):
    sent = {}

    def send_request(eval_id, phase, theta, **kw):
        sent[eval_id] = list(theta)

    def wait_for_result(eval_id, results_path, failures_path, **kw):
        if stop_after is not None and eval_id > stop_after:
            raise Stop
        time.sleep(MATLAB_SECONDS)
        if eval_id in infeasible_ids:
            # What serve_requests writes: a failures row carrying the identifier.
            # append_timeout_failure always writes 'timeout', so it is restamped.
            mi.append_timeout_failure(phase_key, eval_id, "20261002_000000",
                                      "step 3 of case 1 used 70.2 s against a 66.0 s deadline",
                                      sent[eval_id])
            text = Path(failures_path).read_text().replace(",timeout,", f",{RT_ID},")
            Path(failures_path).write_text(text)
            raise mi.EvaluationFailed(eval_id, {"eval_id": eval_id, "identifier": RT_ID,
                                                "message": "deadline", "timestamp": ""})
        row = measured_row(eval_id, phase_label, sent[eval_id], 1000.0 + 10 * eval_id)
        mi.append_imputed_result(phase_key, row)
        return [r for r in mi.read_results(results_path) if r["eval_id"] == eval_id][0]

    driver.send_request = send_request
    driver.wait_for_result = wait_for_result
    driver.wait_for_matlab_ready = lambda *a, **k: None


cfg = run_config.RunConfig(case="baseline")

# The manifest records the tolerance
check(getattr(cfg, "realtime_tolerance", None) == 1.10, "run_config declares realtime_tolerance = 1.10")
check(cfg.to_dict().get("realtime_tolerance") == 1.10, "the manifest carries realtime_tolerance")

# Design phase: point 3 misses the deadline
make_fakes("init", "DOE", infeasible_ids={3})
try:
    driver.run_initialization(cfg)
    init_ok = True
except Exception as exc:  # noqa: BLE001
    init_ok = False
    print(f"run_initialization raised: {exc!r}")
check(init_ok, "run_initialization completes past a real-time abort")
rows = mi.read_results(mi.results_file("init"))
r3 = [r for r in rows if r["eval_id"] == 3]
check(len(rows) == cfg.n_init, f"design has {len(rows)} rows, expected {cfg.n_init}")
check(len(r3) == 1, "the infeasible design point is in results.csv")
if r3:
    worst = max(r["SSE"] for r in rows if r["eval_id"] < 3)
    check(abs(r3[0]["SSE"] - 1.1 * worst) < 1e-6,
          f"imputed SSE {r3[0]['SSE']} = 1.1 * worst so far {worst}")
    check(MATLAB_SECONDS * 0.9 <= r3[0]["runtime_s"] < 60.0,
          f"imputed runtime_s is the elapsed time ({r3[0]['runtime_s']:.2f} s), not eval_timeout_s")
failure_rows = [f for f in mi.read_failures(mi.failures_file("init")) if f["eval_id"] == 3]
check([f["identifier"] for f in failure_rows] == [RT_ID],
      f"one failures row, MATLAB's, for the point: {[f['identifier'] for f in failure_rows]}")

# Optimization phase: the first BO proposal misses the deadline
make_fakes("bo", "OPT", infeasible_ids={cfg.n_init + 1}, stop_after=cfg.n_init + 1)
try:
    driver.run_bo(cfg)
except Stop:
    pass
bo_rows = mi.read_results(mi.results_file("bo"))
check([r["eval_id"] for r in bo_rows] == [cfg.n_init + 1],
      f"BO history holds the imputed row: {[r['eval_id'] for r in bo_rows]}")
if bo_rows:
    worst = max(r["SSE"] for r in rows)
    check(abs(bo_rows[0]["SSE"] - 1.1 * worst) < 1e-6, "imputed BO SSE is 1.1 * worst observed")
    check(bo_rows[0]["runtime_s"] < 60.0,
          f"imputed BO runtime_s is elapsed ({bo_rows[0]['runtime_s']:.2f} s)")
ledger_path = WORK / "registry" / "evaluations_bo.jsonl"
ledger = ledger_path.read_text() if ledger_path.exists() else ""
check('"imputed_realtime": true' in ledger, "the ledger marks the point imputed_realtime")

shutil.rmtree(WORK, ignore_errors=True)
finish()
