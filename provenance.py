"""Durable record of everything a run did.

Three artefacts, all under results/registry:

  manifest.json          written once, before the first evaluation. Holds the
                         declared configuration, the environment, the git state
                         and the seeds. On resume it is compared against the
                         current configuration and any difference stops the run.

  evaluations.jsonl      one line per evaluation, appended and flushed to disk
                         before the next request is sent. Holds the proposal
                         (acquisition value, reference point, fitted GP
                         hyperparameters, surrogate vintage in force) and the
                         result read back from MATLAB.

  vintages/vintage_NN.json   one file per surrogate fit: shape parameters,
                         selected lambda, the whole cross-validation loss grid,
                         the optimiser exit status, and the list of output files
                         that informed the fit.

Design rules, both consequences of expecting a crash at any point:

  Append, never rewrite. A line already on disk is never revised, so a partial
  write can only ever truncate the tail. The reader drops a trailing line that
  does not parse and the run resumes from the last complete record.

  Write the durable record after the fact it describes. The evaluation line is
  appended once MATLAB has reported the result, so a line in the ledger implies
  a completed evaluation. The reverse gap, a completed evaluation with no line,
  is recoverable from the results CSV and is repaired on resume.
"""

from __future__ import annotations

import json
import os
import platform
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional


REGISTRY_DIRNAME = "registry"


# ---------------------------------------------------------------------------
# Durable primitives
# ---------------------------------------------------------------------------

def atomic_write_json(path: Path, payload: Dict) -> None:
    """Write a JSON file by temporary name and rename.

    A reader sees either the previous content or the new content. The fsync
    before the rename is what makes that true after a power loss rather than
    only after a process kill.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True, default=str)
        fh.flush()
        os.fsync(fh.fileno())
    _replace_with_retry(tmp, path)


def _replace_with_retry(src: Path, dst: Path, attempts: int = 10,
                        pause_s: float = 0.15) -> None:
    """Rename over a destination another process may have open.

    On Windows a rename fails while the target is open for reading, which the
    MATLAB side does whenever it reloads coefficients. The collision window is
    a few milliseconds, so a bounded retry closes it.
    """
    last: Optional[Exception] = None
    for _ in range(attempts):
        try:
            src.replace(dst)
            return
        except OSError as exc:
            last = exc
            time.sleep(pause_s)
    raise RuntimeError(f"could not replace {dst} after {attempts} attempts") from last


def append_jsonl(path: Path, record: Dict) -> None:
    """Append one JSON object as a line, flushed to disk before returning."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    line = json.dumps(record, sort_keys=True, default=str)
    with path.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")
        fh.flush()
        os.fsync(fh.fileno())


def read_jsonl(path: Path) -> List[Dict]:
    """Read a JSONL file, discarding a truncated final line.

    A crash during the append leaves a partial last line. Dropping it is
    correct: the record it described is re-derived from the results CSV on
    resume, and keeping it would poison every later parse.
    """
    path = Path(path)
    if not path.exists():
        return []

    records: List[Dict] = []
    with path.open("r", encoding="utf-8") as fh:
        lines = fh.readlines()

    for i, line in enumerate(lines):
        line = line.strip()
        if not line:
            continue
        try:
            records.append(json.loads(line))
        except json.JSONDecodeError:
            if i == len(lines) - 1:
                print(f"[registry] dropping a truncated final line in {path.name}")
                break
            raise
    return records


# ---------------------------------------------------------------------------
# Environment capture
# ---------------------------------------------------------------------------

def _git(*args: str, cwd: Path) -> Optional[str]:
    try:
        out = subprocess.run(
            ["git", *args], cwd=str(cwd), capture_output=True, text=True, timeout=15)
        if out.returncode != 0:
            return None
        return out.stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return None


_ENV_CACHE: Dict[str, tuple] = {}


def capture_environment(root: Path) -> tuple:
    """Record what the run was executed with.

    The dirty flag and the diff of tracked source files matter more than the
    commit hash: a run started from an edited working tree is not reproducible
    from the hash alone, and recording the diff is what makes it reproducible
    anyway.

    Cached per process. The git calls cost seconds on a network-backed working
    tree, and nothing they report can change while one run is executing.
    """
    key = str(Path(root).resolve())
    if key in _ENV_CACHE:
        return _ENV_CACHE[key]

    env: Dict = {
        "python": sys.version,
        "executable": sys.executable,
        "platform": platform.platform(),
        "processor": platform.processor(),
        "machine": platform.machine(),
        "cwd": str(Path.cwd()),
        "started_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "started_at_epoch": time.time(),
    }

    packages: Dict[str, str] = {}
    for name in ("torch", "botorch", "gpytorch", "numpy", "scipy"):
        try:
            module = __import__(name)
            packages[name] = getattr(module, "__version__", "unknown")
        except ImportError:
            packages[name] = "not installed"
    env["packages"] = packages

    try:
        import torch
        env["torch_device"] = "cuda" if torch.cuda.is_available() else "cpu"
        env["torch_num_threads"] = torch.get_num_threads()
        if torch.cuda.is_available():
            env["cuda_device_name"] = torch.cuda.get_device_name(0)
            env["cuda_version"] = torch.version.cuda
    except ImportError:
        pass

    # The diff is scoped to source files. The repository also tracks simulation
    # outputs, and diffing those costs seconds per start while saying nothing
    # about what the run will do. Their state is recorded by the results files
    # themselves.
    source_globs = ["*.py", "*.m", "*.json", "*.toml", "*.cfg"]
    diff = _git("diff", "HEAD", "--", *source_globs, cwd=root)
    changed = _git("diff", "HEAD", "--name-only", "--", *source_globs, cwd=root)

    env["git"] = {
        "commit": _git("rev-parse", "HEAD", cwd=root),
        "branch": _git("rev-parse", "--abbrev-ref", "HEAD", cwd=root),
        "describe": _git("describe", "--always", "--dirty", cwd=root),
        "dirty": bool(diff),
        "changed_source_files": changed.splitlines() if changed else [],
        "diff_bytes": len(diff) if diff else 0,
        "diff_scope": source_globs,
    }

    _ENV_CACHE[key] = (env, diff)
    return env, diff


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

@dataclass
class Registry:
    """Paths and operations of one run's durable record."""

    root: Path
    phase: str

    def __post_init__(self) -> None:
        self.dir = Path(self.root) / REGISTRY_DIRNAME
        self.dir.mkdir(parents=True, exist_ok=True)
        self.manifest_path = self.dir / f"manifest_{self.phase}.json"
        self.evaluations_path = self.dir / f"evaluations_{self.phase}.jsonl"
        self.vintages_dir = self.dir / "vintages"
        self.diff_path = self.dir / f"working_tree_{self.phase}.diff"

    # -- manifest ----------------------------------------------------------

    def write_manifest(self, config: Dict, project_root: Path, extra: Dict | None = None) -> Dict:
        """Write the manifest, or verify the existing one against this run.

        A resumed run must not change anything the manifest declares. Comparing
        the configuration block verbatim is what turns a silent divergence, an
        edited budget or a different case, into an error at start rather than an
        inconsistency discovered in the results months later.
        """
        env, diff = capture_environment(project_root)
        manifest = {
            "phase": self.phase,
            "config": config,
            "environment": env,
            "registry_version": 2,
        }
        if extra:
            manifest["context"] = extra

        if self.manifest_path.exists():
            previous = json.loads(self.manifest_path.read_text(encoding="utf-8"))
            self._assert_config_matches(previous.get("config", {}), config)
            resumes = previous.setdefault("resumed", [])
            resumes.append({
                "at": env["started_at"],
                "git_commit": env["git"]["commit"],
                "git_dirty": env["git"]["dirty"],
            })
            atomic_write_json(self.manifest_path, previous)
            print(f"[registry] resuming under manifest {self.manifest_path.name} "
                  f"({len(resumes)} restart(s) recorded)")
            return previous

        atomic_write_json(self.manifest_path, manifest)
        if diff:
            self.diff_path.write_text(diff, encoding="utf-8")
            print(f"[registry] working tree was dirty; diff saved to {self.diff_path.name}")
        print(f"[registry] manifest written to {self.manifest_path}")
        return manifest

    @staticmethod
    def _assert_config_matches(previous: Dict, current: Dict) -> None:
        # Compare both sides after a round trip through JSON. The stored copy has
        # already been through it, so a tuple in the live configuration would
        # otherwise differ from the list it was written as, and an identical run
        # would be refused on resume.
        previous = json.loads(json.dumps(previous, sort_keys=True, default=str))
        current = json.loads(json.dumps(current, sort_keys=True, default=str))

        differences = []
        for key in sorted(set(previous) | set(current)):
            was, now = previous.get(key, "<absent>"), current.get(key, "<absent>")
            if was != now:
                differences.append(f"  {key}: manifest has {was!r}, run has {now!r}")
        if differences:
            raise RuntimeError(
                "the configuration differs from the manifest this run started under:\n"
                + "\n".join(differences)
                + "\n\nResuming under a changed configuration would mix two runs in one "
                  "results file. Restore the setting, or start in a fresh results "
                  "directory."
            )

    # -- evaluations -------------------------------------------------------

    def append_evaluation(self, record: Dict) -> None:
        append_jsonl(self.evaluations_path, record)

    def load_evaluations(self) -> List[Dict]:
        return read_jsonl(self.evaluations_path)

    def last_eval_id(self) -> int:
        records = self.load_evaluations()
        if not records:
            return 0
        return max(int(r.get("eval_id", 0)) for r in records)

    def reconcile(self, rows: List[Dict]) -> int:
        """Fold results rows the ledger is missing back into it.

        The gap this repairs is a crash between MATLAB appending its row and the
        driver appending its own. The evaluation happened and its result is on
        disk, so it belongs in the history; what is lost is the proposal
        metadata, which the recovered entry marks as absent rather than
        inventing.

        Idempotent, so it can run on every start. Returns the number of entries
        recovered.
        """
        known = {int(r.get("eval_id", 0)) for r in self.load_evaluations()}
        recovered = 0
        for row in rows:
            if int(row["eval_id"]) in known:
                continue
            self.append_evaluation({
                "eval_id": int(row["eval_id"]),
                "recovered": True,
                "note": "reconstructed from the results CSV; the proposal metadata "
                        "was lost to an interruption between the MATLAB write and "
                        "the driver write",
                "result": row,
                "recorded_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            })
            known.add(int(row["eval_id"]))
            recovered += 1
        if recovered:
            print(f"[registry] recovered {recovered} evaluation(s) from the results CSV")
        return recovered

    # -- vintages ----------------------------------------------------------

    def vintage_path(self, vintage: int) -> Path:
        return self.vintages_dir / f"vintage_{vintage:02d}.json"

    def vintage_exists(self, vintage: int) -> bool:
        return self.vintage_path(vintage).exists()

    def load_vintage(self, vintage: int) -> Optional[Dict]:
        path = self.vintage_path(vintage)
        if not path.exists():
            return None
        return json.loads(path.read_text(encoding="utf-8"))

    def list_vintages(self) -> List[int]:
        if not self.vintages_dir.exists():
            return []
        out = []
        for p in sorted(self.vintages_dir.glob("vintage_*.json")):
            try:
                out.append(int(p.stem.split("_")[1]))
            except (IndexError, ValueError):
                continue
        return out


def summarise_gp(model, name: str) -> Dict:
    """Read the fitted hyperparameters of one GP into a plain dictionary.

    Recorded per proposal, because they are what turned the data into the
    surface the acquisition was maximised over. Without them a proposal cannot
    be reproduced from the data alone: the marginal likelihood has local optima
    and the fit is warm-started from whatever the previous iteration reached.
    """
    import torch

    out: Dict = {"name": name}
    try:
        with torch.no_grad():
            covar = model.covar_module
            if hasattr(covar, "base_kernel") and hasattr(covar.base_kernel, "lengthscale"):
                out["lengthscale"] = covar.base_kernel.lengthscale.detach().cpu().flatten().tolist()
            elif hasattr(covar, "lengthscale"):
                out["lengthscale"] = covar.lengthscale.detach().cpu().flatten().tolist()
            if hasattr(covar, "outputscale"):
                out["outputscale"] = covar.outputscale.detach().cpu().flatten().tolist()
            if hasattr(model.likelihood, "noise"):
                out["noise"] = model.likelihood.noise.detach().cpu().flatten().tolist()
            if hasattr(model, "mean_module") and hasattr(model.mean_module, "constant"):
                out["mean_constant"] = model.mean_module.constant.detach().cpu().flatten().tolist()
    except Exception as exc:                      # diagnostics must not stop a run
        out["error"] = f"{type(exc).__name__}: {exc}"
    return out
