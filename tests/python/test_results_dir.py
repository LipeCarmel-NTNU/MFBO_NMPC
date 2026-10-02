"""Per-campaign results folders derived from the case and the Sobol seed.

    python tests/python/test_results_dir.py

Each check starts a fresh interpreter, because the results folder is fixed when
pipeline.matlab_interface is imported.
"""
import os
import re
import subprocess
import sys

from _check import REPO, check, find_bash, finish

PY = sys.executable


def run(args, env_extra=None, unset=True):
    env = dict(os.environ)
    if unset:
        env.pop("MFBO_RESULTS_DIR", None)
    env.update(env_extra or {})
    return subprocess.run([PY, *args], cwd=REPO, env=env, capture_output=True, text=True)


# Names
r = run(["-c", "from run_config import results_dir_for as f; "
               "print(f('baseline')); print(f('case2')); print(f('case1'))"])
check(r.stdout.split() == ["results/SF_baseline_s123", "results/MF_case2_s123",
                           "results/MF_case1_s123"],
      f"results_dir_for names the arm, case and seed: {r.stdout or r.stderr}")

# Derived at import from --case, and exported for the MATLAB child
probe = ("import os, pipeline.matlab_interface as mi; "
         "print(mi.RESULTS_DIR.relative_to(mi.BASE_DIR).as_posix()); "
         "print(os.environ.get('MFBO_RESULTS_DIR'))")
r = run(["-c", probe, "--case", "baseline"])
check(r.stdout.split() == ["results/SF_baseline_s123"] * 2,
      f"--case baseline is derived and exported: {r.stdout or r.stderr}")
r = run(["-c", probe, "--case=case2"])
check(r.stdout.split()[:1] == ["results/MF_case2_s123"], f"--case=case2 form: {r.stdout or r.stderr}")

# The variable still overrides
r = run(["-c", probe, "--case", "baseline"], {"MFBO_RESULTS_DIR": "results/elsewhere"}, unset=False)
check(r.stdout.split()[:1] == ["results/elsewhere"], f"MFBO_RESULTS_DIR overrides: {r.stdout or r.stderr}")

# One case per run_supervised
r = run(["-c", "import run_supervised as s; s.parse(['--case', 'case1', 'case2'])", "--case", "case1"])
check(r.returncode == 2 and "unrecognized arguments: case2" in r.stderr,
      f"run_supervised refuses two cases: {r.stderr[-200:]}")
r = run(["-c", "import run_supervised as s; print(s.parse(['--case', 'baseline']).case)", "--case", "baseline"])
check(r.stdout.strip() == "baseline", f"run_supervised takes one case: {r.stdout or r.stderr[-300:]}")

# watch_matlab_log resolves the campaign log from --case
r = run(["watch_matlab_log.py", "--case", "baseline", "--no-follow"])
check("SF_baseline_s123" in (r.stdout + r.stderr), f"watch_matlab_log --case: {(r.stdout + r.stderr)[-200:]}")

# No pilot default left in the code that runs a campaign
code = list(REPO.glob("*.py")) + list(REPO.glob("pipeline/*.py")) + \
       list(REPO.glob("*.m")) + [REPO / "run_idun.slurm"]
hits = [p.name for p in code
        if p.name != "main_rdu_damping.m" and "case2_v3" in p.read_text(errors="replace")]
check(not hits, f"no case2_v3 default in {hits}")

# The archive machinery is gone
check(not (REPO / "pipeline/case_archive.py").exists(), "case_archive.py removed")
check("case_archive" not in (REPO / "pipeline/version.py").read_text(), "version.SOURCES drops case_archive")

# MATLAB servers stop when the variable is missing
check('"MFBO:resultsDir"' in (REPO / "dependencies/io/server_config.m").read_text(),
      "server_config.m raises MFBO:resultsDir")
for m in ("main_initialization.m", "main_BO.m"):
    check("server_config(" in (REPO / m).read_text(), f"{m} takes its folder from server_config")

# Slurm: one case, pre-flight on the derived tree
slurm = (REPO / "run_idun.slurm").read_text()
check("run_idun.slurm all" not in slurm, "slurm drops the all mode")
check('"$MFBO_RESULTS_DIR/init/results.csv"' in slurm, "slurm bo pre-flight checks init/results.csv")
check(not re.search(r"(?<![\w$/])results/(logs|surrogate)", slurm),
      "slurm has no hard-coded results/logs or results/surrogate")
bash = find_bash()
r = subprocess.run([bash, "-n", str(REPO / "run_idun.slurm")], capture_output=True, text=True) if bash else None
check(r is not None and r.returncode == 0, f"slurm syntax: {r.stderr if r else 'no bash found'}")

finish()
