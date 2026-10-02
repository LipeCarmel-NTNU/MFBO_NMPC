"""Pass/fail bookkeeping shared by the test scripts in this folder.

Each test script is a plain program: it prints one PASS or FAIL line per check
and exits non-zero when any check failed. run_all.py runs every script.
"""
import shutil
import sys
from pathlib import Path
from typing import List, Optional

REPO = Path(__file__).resolve().parents[2]
FAILS: List[str] = []


def check(cond: bool, msg: str) -> None:
    print(("PASS " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


def finish() -> None:
    print(f"\n{len(FAILS)} failure(s)")
    sys.exit(1 if FAILS else 0)


def find_bash() -> Optional[str]:
    """A POSIX bash. On Windows, the one beside git, not the WSL launcher."""
    git = shutil.which("git")
    if git:
        root = Path(git).resolve().parents[1]
        for candidate in (root / "usr" / "bin" / "bash.exe", root / "bin" / "bash.exe"):
            if candidate.exists():
                return str(candidate)
    found = shutil.which("bash")
    if found and "system32" not in found.lower():
        return found
    return None
