"""Settings that several files carry must agree.

    python tests/python/test_settings.py

The worker count lives in server_config.m and in run_idun.slurm, and the lock
staleness in run_config.py and server_config.m. A disagreement between the
copies changes no line of code and fails in production, so it is checked here.
"""
import re

from _check import REPO, check, finish

server = (REPO / "dependencies/io/server_config.m").read_text()
slurm = (REPO / "run_idun.slurm").read_text()


def matlab_value(name):
    m = re.search(rf"^\s*cfg_run\.{name}\s*=\s*([^;]+);", server, re.M)
    return eval(m.group(1)) if m else None  # the values are plain arithmetic


# Workers: 8 MATLAB workers, and the reservation one core larger
workers = matlab_value("NumWorkers")
slurm_workers = re.search(r"^MATLAB_WORKERS=(\d+)", slurm, re.M)
cpus = re.search(r"^#SBATCH --cpus-per-task=(\d+)$", slurm, re.M)
check(workers == 8, f"server_config.m opens 8 workers ({workers})")
check(slurm_workers and int(slurm_workers.group(1)) == workers,
      "run_idun.slurm MATLAB_WORKERS matches server_config.m")
check(cpus and int(cpus.group(1)) == workers + 1, "--cpus-per-task is the worker count plus one")
for m in ("main_initialization.m", "main_BO.m"):
    src = (REPO / m).read_text()
    check("configure_pool(true, cfg_run.NumWorkers)" in src, f"{m} opens the server_config pool")
    check(not re.search(r"^NumWorkers\s*=", src, re.M), f"{m} sets no worker count of its own")

# Lock staleness: one value, and at least the evaluation timeout
import sys  # noqa: E402
sys.path.insert(0, str(REPO))
sys.argv = [sys.argv[0], "--case", "case2"]
import inspect  # noqa: E402
import run_config  # noqa: E402
from pipeline import matlab_interface as mi  # noqa: E402

cfg = run_config.RunConfig()
check(cfg.lock_stale_s >= cfg.eval_timeout_s, "lock_stale_s is not below eval_timeout_s")
check(inspect.signature(mi.send_request).parameters["lock_stale_s"].default == cfg.lock_stale_s,
      "send_request defaults to RunConfig.lock_stale_s")
check(inspect.signature(mi.wait_for_lock).parameters["stale_s"].default == cfg.lock_stale_s,
      "wait_for_lock defaults to RunConfig.lock_stale_s")
check(matlab_value("lock_stale_s") == cfg.lock_stale_s, "server_config.m mirrors lock_stale_s")
check(matlab_value("max_consecutive_failures") == cfg.max_consecutive_failures,
      "server_config.m mirrors max_consecutive_failures")

# The real-time tolerance: nmpc_base.m is live, run_config.py records it
base_src = (REPO / "dependencies/simulation/nmpc_base.m").read_text()
tol = re.search(r"^\s*base\.rt_tolerance\s*=\s*([\d.]+);", base_src, re.M)
check(tol and float(tol.group(1)) == cfg.realtime_tolerance, "nmpc_base.m and run_config.py agree on the tolerance")

finish()
