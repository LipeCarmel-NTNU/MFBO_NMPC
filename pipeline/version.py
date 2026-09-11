"""What code is about to run, printed before it runs.

A run once carried a run_config.py that declared doe_prefix_z while the
driver.py beside it had no idea what that setting was. The manifest recorded the
setting, the config was honoured everywhere it was read, and the feature it named
never executed. Nothing in the log said so. Twenty hours of cluster time went
into answering a question the code could not ask.

This module makes that failure visible in the first three lines of every run. It
reports two things.

A digest of each source file that decides behaviour, and one combined digest over
all of them. The digest is taken over the file with newlines normalised, so a
Windows checkout and a Linux checkout of the same commit agree. Compare the
combined digest here against the one on the machine you trust; if they match,
the two are running the same pipeline, whatever git says about branches and
dirty trees.

A capability line, read out of the source rather than out of a constant someone
has to remember to bump. A feature appears as "on" only when the function that
implements it is defined and the call site that uses it is present. That is the
check that would have caught the stale driver: run_config said doe_prefix on,
the driver had no prefix call site, and the two would have disagreed in print.

The features are read by parsing the files, not by importing them, so this
module needs nothing beyond the standard library. It therefore answers on a
login node with no environment activated, which is when you want to ask.

Standalone, to check a machine before submitting anything to it:

    python -m pipeline.version
"""

from __future__ import annotations

import ast
import hashlib
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

BASE_DIR = Path(__file__).resolve().parents[1]

# The files whose content changes what a run does. A file absent from this list
# is a file whose edits this banner will not notice, so add one here when it
# starts carrying behaviour.
SOURCES: Tuple[str, ...] = (
    "run_config.py",
    "run_pipeline.py",
    "pipeline/driver.py",
    "pipeline/fidelity_rows.py",
    "pipeline/phi_surrogate.py",
    "pipeline/matlab_interface.py",
    "pipeline/matlab_supervisor.py",
    "pipeline/provenance.py",
    "pipeline/case_archive.py",
    "pipeline/console_log.py",
)

_DIGEST_CHARS = 8


def file_digest(path: Path) -> Optional[str]:
    """Short content digest, newline-normalised so checkouts compare equal."""
    try:
        raw = path.read_bytes()
    except OSError:
        return None
    normalised = raw.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(normalised).hexdigest()[:_DIGEST_CHARS]


def module_report(base: Path = BASE_DIR) -> List[Dict]:
    """One record per source file: digest, size and modification time."""
    out: List[Dict] = []
    for rel in SOURCES:
        path = base / rel
        digest = file_digest(path)
        if digest is None:
            out.append({"file": rel, "digest": None, "bytes": None, "mtime": None})
            continue
        stat = path.stat()
        out.append({
            "file": rel,
            "digest": digest,
            "bytes": stat.st_size,
            "mtime": time.strftime("%Y-%m-%d %H:%M", time.localtime(stat.st_mtime)),
        })
    return out


def combined_digest(records: List[Dict]) -> str:
    """One digest over every source file, in the declared order.

    A missing file contributes the literal word "missing", so a pipeline with
    fidelity_rows.py absent cannot collide with one that has it.
    """
    h = hashlib.sha256()
    for rec in records:
        h.update(rec["file"].encode("utf-8"))
        h.update(b"\0")
        h.update((rec["digest"] or "missing").encode("utf-8"))
        h.update(b"\0")
    return h.hexdigest()[:12]


def _parse(path: Path):
    """Parse a source file, or None when it is absent or does not parse."""
    try:
        return ast.parse(path.read_text(encoding="utf-8"))
    except (OSError, SyntaxError):
        return None


def _defines(tree, name: str):
    """The top-level function definition called name, or None."""
    if tree is None:
        return None
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name:
            return node
    return None


def _assigns(tree, name: str) -> bool:
    """Whether the module assigns name at the top level."""
    if tree is None:
        return False
    for node in tree.body:
        targets = []
        if isinstance(node, ast.Assign):
            targets = node.targets
        elif isinstance(node, ast.AnnAssign):
            targets = [node.target]
        if any(isinstance(t, ast.Name) and t.id == name for t in targets):
            return True
    return False


def _mentions(node, name: str) -> bool:
    """Whether name appears anywhere inside a function body."""
    if node is None:
        return False
    return any(isinstance(n, ast.Name) and n.id == name for n in ast.walk(node))


def capabilities(base: Path = BASE_DIR) -> Dict[str, str]:
    """What the source on this machine can actually do.

    Read by parsing, so the answer does not depend on torch, botorch or any
    environment being importable, and a file that will not parse is reported as
    such rather than taken for a file without the feature.
    """
    driver_tree = _parse(base / "pipeline/driver.py")
    if driver_tree is None:
        return {"driver_source": "UNREADABLE (missing or does not parse)"}

    caps: Dict[str, str] = {}

    # Design-prefix rows need both halves: the reader, and the call in run_bo
    # that puts its rows into the history. The run that prompted this module had
    # neither, while its configuration named the fidelities.
    reader = _defines(_parse(base / "pipeline/fidelity_rows.py"), "doe_prefix_rows")
    run_bo = _defines(driver_tree, "run_bo")
    caps["doe_prefix"] = _both(reader is not None,
                               _mentions(run_bo, "doe_prefix_rows"),
                               "fidelity_rows.doe_prefix_rows", "call in run_bo")

    # Re-measuring the whole history under the vintage in force is the second
    # argument of history_tensors. Without it the design rows enter the
    # objective GP unextrapolated while the optimisation rows are scaled.
    ht = _defines(driver_tree, "history_tensors")
    caps["rescale_past_rows"] = "on" if (
        ht is not None and len(ht.args.args) >= 2) else "off"

    caps["phi_floor"] = "on" if _assigns(driver_tree, "PHI_FLOOR") else "off"
    caps["timeout_impute"] = "on" if _defines(driver_tree, "_impute_timeout") else "off"
    return caps


def _both(a: bool, b: bool, name_a: str, name_b: str) -> str:
    if a and b:
        return "on"
    if not a and not b:
        return "off"
    return f"BROKEN (no {name_b if a else name_a})"


def git_line(base: Path = BASE_DIR) -> str:
    """Branch, description and dirty state, or a note that git is unavailable."""
    try:
        from pipeline.provenance import _git
    except Exception:                                         # noqa: BLE001
        return "unavailable"
    branch = _git("rev-parse", "--abbrev-ref", "HEAD", cwd=base) or "?"
    describe = _git("describe", "--always", "--dirty", cwd=base) or "?"
    return f"{branch} @ {describe}"


def report(cfg=None, base: Path = BASE_DIR) -> Dict:
    """Everything the banner prints, as data."""
    records = module_report(base)
    out: Dict = {
        "pipeline_digest": combined_digest(records),
        "modules": records,
        "capabilities": capabilities(base),
        "git": git_line(base),
    }
    if cfg is not None:
        out["settings"] = {
            "doe_prefix_z": list(getattr(cfg, "doe_prefix_z", ()) or ()),
            "refit_every": getattr(cfg, "refit_every", None),
            "is_baseline": getattr(cfg, "is_baseline", None),
        }
    return out


def print_banner(cfg=None, base: Path = BASE_DIR) -> Dict:
    """Print the banner and return the same content as data."""
    rep = report(cfg, base)

    print(f"[ver] pipeline {rep['pipeline_digest']}   git {rep['git']}")
    for rec in rep["modules"]:
        if rec["digest"] is None:
            print(f"[ver]   {rec['file']:<32} MISSING")
        else:
            print(f"[ver]   {rec['file']:<32} {rec['digest']}  "
                  f"{rec['bytes']:>7} B  {rec['mtime']}")

    caps = rep["capabilities"]
    print("[ver] features: " + "  ".join(f"{k}={v}" for k, v in caps.items()))

    settings = rep.get("settings")
    if settings is not None:
        zs = settings["doe_prefix_z"]
        print(f"[ver] settings: doe_prefix_z="
              f"{'(' + ', '.join(f'{z:g}' for z in zs) + ')' if zs else 'none'}  "
              f"refit_every={settings['refit_every']}  "
              f"is_baseline={settings['is_baseline']}")

        # The exact contradiction that went unnoticed: the configuration asks
        # for prefix rows and the code cannot produce them. Say so loudly, here,
        # rather than leaving it to be inferred from a missing log line later.
        if zs and not settings["is_baseline"] and caps.get("doe_prefix") != "on":
            print("[ver] *** doe_prefix_z is set but this pipeline cannot use it "
                  f"(doe_prefix={caps.get('doe_prefix')}). The design prefix rows "
                  "will NOT be added. Update pipeline/ before running. ***")

    return rep


if __name__ == "__main__":
    print_banner()
