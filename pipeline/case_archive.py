"""Move a finished case out of results/ so that the next case can start there.

Both cases write the same paths. main_initialization.m and main_BO.m name
results/init and results, and the driver refuses to resume under a manifest that
declares a different case, so two cases cannot share one results tree. RERUN.md
states the manual procedure: move results/ aside and repeat with the other case.
This module is that step, performed between cases by run_supervised.py.

The move is a rename, so nothing is copied and nothing is deleted. The
destination carries the case name and a timestamp and is never reused, which
means a mistake here leaves the earlier run recoverable by moving the directory
back.
"""

from __future__ import annotations

import json
import shutil
import time
from pathlib import Path
from typing import List, Optional

BASE_DIR = Path(__file__).resolve().parents[1]
RESULTS_DIR = BASE_DIR / "results"
ARCHIVE_ROOT = BASE_DIR / "results_archive"
SIM_LOG = BASE_DIR / "SIMULATIONS_LOG.txt"

# The MATLAB scripts write these two at the top of the working directory rather
# than under results/, so an archive of results/ leaves them behind.
CARRIED_FILES = (SIM_LOG,)

# The exchange has to be empty before another case starts. A request left from
# the previous case names an eval_id that the new driver has not sent, and a
# freshly started server has no record of having answered it, so it would be
# served and appended to the new results file.
THETA_FILE = BASE_DIR / "inbox" / "theta.txt"
LOCK_FILE = BASE_DIR / "matlab.lock"


def declared_case(results_dir: Path = RESULTS_DIR) -> Optional[str]:
    """Return the case that the results tree belongs to, or None.

    The manifest is the only record that names the case, and the driver compares
    against it on resume. Reading it here is what lets the sequence tell a
    resumable tree from one that has to be moved.
    """
    registry = Path(results_dir) / "registry"
    for name in ("manifest_init.json", "manifest_bo.json"):
        path = registry / name
        if not path.exists():
            continue
        try:
            manifest = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        case = manifest.get("config", {}).get("case")
        if isinstance(case, str) and case:
            return case
    return None


def is_empty(results_dir: Path = RESULTS_DIR) -> bool:
    path = Path(results_dir)
    if not path.exists():
        return True
    try:
        return not any(path.iterdir())
    except OSError:
        return False


def archive(results_dir: Path = RESULTS_DIR, *,
            archive_root: Path = ARCHIVE_ROOT,
            case_hint: Optional[str] = None,
            carried: tuple = CARRIED_FILES) -> Optional[Path]:
    """Move the results tree into archive_root and leave an empty one behind.

    Returns the destination, or None when there was nothing to move.
    """
    results_dir = Path(results_dir)
    if is_empty(results_dir):
        results_dir.mkdir(parents=True, exist_ok=True)
        return None

    case = declared_case(results_dir) or case_hint or "unknown_case"
    archive_root = Path(archive_root)
    archive_root.mkdir(parents=True, exist_ok=True)

    dest = archive_root / f"{case}_{time.strftime('%Y%m%d_%H%M%S')}"
    if dest.exists():
        # Two archives inside one second. The suffix keeps the earlier one.
        for n in range(2, 100):
            candidate = archive_root / f"{case}_{time.strftime('%Y%m%d_%H%M%S')}_{n}"
            if not candidate.exists():
                dest = candidate
                break
        else:
            raise RuntimeError(f"could not find an unused archive name under {archive_root}")

    shutil.move(str(results_dir), str(dest))
    results_dir.mkdir(parents=True, exist_ok=True)

    for path in carried:
        path = Path(path)
        if path.exists():
            shutil.move(str(path), str(dest / path.name))

    print(f"[archive] {results_dir.name}/ held {case} and was moved to {dest}")
    return dest


def clear_exchange() -> List[Path]:
    """Remove a request and a lock left by the previous case."""
    removed = []
    for path in (THETA_FILE, LOCK_FILE):
        if path.exists():
            try:
                path.unlink()
                removed.append(path)
            except OSError:
                pass
    if removed:
        print("[archive] cleared " + ", ".join(p.name for p in removed))
    return removed


def prepare_for(case: str, *, results_dir: Path = RESULTS_DIR,
                archive_root: Path = ARCHIVE_ROOT) -> Optional[Path]:
    """Make results/ ready to hold this case, and report what was moved.

    A tree that already declares this case is left alone, because that is what a
    resumed run needs. A tree that declares another case is archived. A tree with
    no manifest is treated as belonging to this case: it holds no evaluation that
    the driver would refuse to resume.
    """
    results_dir = Path(results_dir)
    existing = declared_case(results_dir)
    moved = None
    if existing is not None and existing != case:
        moved = archive(results_dir, archive_root=archive_root, case_hint=existing)
    elif existing == case:
        print(f"[archive] results/ already holds {case}; it will be resumed")
    clear_exchange()
    results_dir.mkdir(parents=True, exist_ok=True)
    return moved
