"""Run every Python test script in this folder.

    python tests/python/run_all.py

Each test_*.py runs in its own interpreter, because the results folder is fixed
when pipeline.matlab_interface is imported. Use the project's environment, the
one that imports torch and botorch. The MATLAB tests are in tests/matlab.
"""
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main() -> int:
    failed = []
    for script in sorted(HERE.glob("test_*.py")):
        print(f"\n== {script.name}", flush=True)
        r = subprocess.run([sys.executable, str(script)], cwd=HERE)
        if r.returncode != 0:
            failed.append(script.name)
    print(f"\n{len(failed)} script(s) failed" + (f": {', '.join(failed)}" if failed else ""))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
