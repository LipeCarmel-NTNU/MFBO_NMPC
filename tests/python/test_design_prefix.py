"""The MF arm stops at start when the design trends files are not all there.

    python tests/python/test_design_prefix.py

The MF clone receives the design from the SF clone as init/results.csv plus
every init/out_<ts>.mat. The .mat files feed the design prefix rows and the
vintage-0 fit, and a missing or unreadable one used to drop out without an
error. The files here are real .mat files in the layout simulate_nmpc saves,
read by the real doe_prefix_rows.
"""
import os
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import scipy.io

from _check import REPO, check, finish

WORK = Path(tempfile.mkdtemp(prefix="prefix_"))
os.environ["MFBO_RESULTS_DIR"] = str(WORK)
sys.path.insert(0, str(REPO))
sys.argv = [sys.argv[0], "--case", "case2"]

import run_config  # noqa: E402
from pipeline import driver  # noqa: E402
from pipeline import matlab_interface as mi  # noqa: E402

N = 601
INIT = WORK / "init"
INIT.mkdir(parents=True)


def design_row(eval_id, imputed=False):
    theta = [1.0, 3.0, 2.0, 0.0, 0.5, -0.5, -1000.0, -1000.0, -1000.0, 0.1, 0.2, 0.3]
    wall = 0.0 if imputed else 0.4      # what _impute_dominated writes, and a measured row
    return {"eval_id": eval_id, "timestamp": f"20991231_{eval_id:06d}", "phase": "DOE",
            "phi_vintage": float("nan"), "z": 1.0, "SSE_measured": 1000.0 + eval_id,
            "SSdU_measured": 1e-3, "phi_SSE": 1.0, "phi_SSdU": 1.0, "SSE": 1000.0 + eval_id,
            "SSdU": 1e-3, "J": 1010.0 + eval_id, "runtime_s": 100.0, "n_flag_not_one": 0.0,
            "phi_floored": 0.0, "wall_total_s": 120.0, "wall_cases_s": 110.0,
            "wall_phi_s": 0.0, "wall_build_s": wall, "wall_save_s": wall, "theta": theta}


def write_trends(row):
    case = {"partial_SSE": np.full(N, 1.0), "partial_SSdU": np.full(N - 1, 1e-6),
            "RUNTIME": np.full(N, 0.1)}
    out = {"N": N, "theta": np.array(row["theta"]),
           "case": np.array([case, case], dtype=object)}
    scipy.io.savemat(str(INIT / f"out_{row['timestamp']}.mat"), {"out": out})


def fresh_design(imputed_ids=()):
    if INIT.exists():
        shutil.rmtree(INIT)
    INIT.mkdir(parents=True)
    rows = [design_row(i, imputed=i in imputed_ids) for i in range(1, 21)]
    for row in rows:
        mi.append_imputed_result("init", row)
        if row["eval_id"] not in imputed_ids:
            write_trends(row)
    return mi.read_results(mi.results_file("init"))


def prefix_rows_or_error(cfg, rows):
    try:
        return driver._design_prefix_rows(cfg, rows), None
    except RuntimeError as exc:
        return None, str(exc)
    except AttributeError as exc:
        return None, f"AttributeError: {exc}"


mf = run_config.RunConfig(case="case2")
n_z = len([z for z in mf.doe_prefix_z if z < 1.0])

# A complete copy
rows = fresh_design()
got, err = prefix_rows_or_error(mf, rows)
check(err is None and len(got) == 20 * n_z, f"complete design gives {20 * n_z} prefix rows ({err or len(got)})")

# A forgotten .mat
(INIT / f"out_{rows[6]['timestamp']}.mat").unlink()
got, err = prefix_rows_or_error(mf, rows)
check(got is None and err and rows[6]["timestamp"] in err and "missing" in err,
      f"a missing trends file stops the run and names it: {err}")

# No .mat at all, the trap the check exists for
for p in INIT.glob("out_*.mat"):
    p.unlink()
got, err = prefix_rows_or_error(mf, rows)
check(got is None and err and "20 design trends file" in err, f"no trends files stops the run: {err}")

# A design point imputed by the real-time guard has no .mat, and is not counted
rows = fresh_design(imputed_ids={3})
got, err = prefix_rows_or_error(mf, rows)
check(err is None and len(got) == 19 * n_z, f"an imputed design row is excluded ({err or len(got)})")

# A file that cannot be read
rows = fresh_design()
scipy.io.savemat(str(INIT / f"out_{rows[0]['timestamp']}.mat"), {"not_out": 1})
got, err = prefix_rows_or_error(mf, rows)
check(got is None and err and "expected" in err, f"an unreadable trends file stops the run: {err}")

# The SF arm reads no trends files
for p in INIT.glob("out_*.mat"):
    p.unlink()
sf = run_config.RunConfig(case="baseline")
got, err = prefix_rows_or_error(sf, rows)
check(err is None and got == [], f"the baseline needs no trends files ({err or got})")

# run_bo goes through the check
src = (REPO / "pipeline" / "driver.py").read_text()
run_bo_src = src[src.index("def run_bo("):src.index("def ", src.index("def run_bo(") + 10)]
check("_design_prefix_rows(cfg, init_rows)" in run_bo_src, "run_bo builds its prefix rows through the check")
from pipeline.version import capabilities  # noqa: E402
check(capabilities()["doe_prefix"] == "on", f"the version banner reports doe_prefix=on ({capabilities()['doe_prefix']})")

shutil.rmtree(WORK, ignore_errors=True)
finish()
