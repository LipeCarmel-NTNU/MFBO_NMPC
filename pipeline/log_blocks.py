"""Mirror the growing MATLAB console log into small dated blocks.

A headless run sends everything the MATLAB command window would have shown to
one appended file, results/<case>/logs/matlab_console.log. Over a campaign that
file reaches tens of megabytes, and reading the last few minutes of it means
opening the whole thing. This module writes a second, block-structured view of
the same bytes beside it:

    logs/matlab_console.log              the raw append-only file, unchanged
    logs/console_blocks/
        block_0000_20260918_142530.log   one block, at most BLOCK_BYTES
        block_0001_20260918_151204.log
        index.csv                        one row per block
        latest.txt                       the name of the block being written
        state.json                       the reader offset, for a clean resume

The raw file is read and never written, so MATLAB's file descriptor is left
alone. That matters on Windows, where renaming or truncating a file another
process holds open is the operation most likely to fail. The cost is a second
copy of a text log, which is megabytes.

A block closes when it reaches BLOCK_BYTES or when it has been open for
BLOCK_SECONDS, whichever comes first. The time limit is what makes the newest
block genuinely recent: without it a quiet period leaves you reading a block
that was opened an hour ago and is still a tenth full.

Blocks split on line boundaries. A line that arrives while a block is full still
lands whole in the next block.

index.csv carries the evaluation ids seen in each block, so finding the block
that holds evaluation 47 is a search of one small file rather than of the log.

Settings come from the environment and not from run_config.py. Every value in
RunConfig reaches the manifest, and the manifest is compared verbatim on resume,
so a logging knob there would refuse a resume after you changed how much output
you wanted to keep.

    MFBO_LOG_BLOCKS=0           turn the mirror off
    MFBO_LOG_BLOCK_BYTES=524288 roll at 512 KiB instead of 1 MiB
    MFBO_LOG_BLOCK_SECONDS=600  roll after 10 minutes instead of 30

Run it standalone against a log the supervisor is not watching:

    python -m pipeline.log_blocks                  # catch up once and exit
    python -m pipeline.log_blocks --follow         # keep mirroring
    python -m pipeline.log_blocks --index          # print index.csv

Only one process may mirror one log. The supervisor runs this on its monitor
thread, so do not also run it standalone against the same file while a
supervised run is going.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional

BLOCKS_DIRNAME = "console_blocks"
STATE_NAME = "state.json"
INDEX_NAME = "index.csv"
LATEST_NAME = "latest.txt"

DEFAULT_BLOCK_BYTES = 1024 * 1024
DEFAULT_BLOCK_SECONDS = 30 * 60.0

INDEX_COLUMNS = ["block", "opened_at", "closed_at", "bytes", "lines",
                 "first_eval", "last_eval", "note"]

# serve_requests.m prints "  eval %d served in %.1f s" and "  eval %d FAILED",
# and both entry points print "  OPT eval %d [phi v%d]" or the DOE equivalent.
# One permissive pattern covers all of them.
EVAL_PATTERN = re.compile(rb"\beval\s+(\d+)\b")


def _env_float(name: str, fallback: float) -> float:
    raw = os.environ.get(name)
    if not raw:
        return fallback
    try:
        return float(raw)
    except ValueError:
        print(f"[blocks] {name}={raw!r} is not a number. Using {fallback}.",
              file=sys.stderr, flush=True)
        return fallback


def blocks_enabled() -> bool:
    return os.environ.get("MFBO_LOG_BLOCKS", "1").strip() not in {"0", "false", "no"}


def block_bytes() -> int:
    return int(_env_float("MFBO_LOG_BLOCK_BYTES", DEFAULT_BLOCK_BYTES))


def block_seconds() -> float:
    return _env_float("MFBO_LOG_BLOCK_SECONDS", DEFAULT_BLOCK_SECONDS)


def blocks_dir_for(source: Path) -> Path:
    """Where the blocks of a given log live."""
    return Path(source).parent / BLOCKS_DIRNAME


def _stamp(when: Optional[float] = None) -> str:
    return time.strftime("%Y%m%d_%H%M%S", time.localtime(when))


def read_index(blocks_dir: Path) -> List[Dict[str, str]]:
    """Every row of index.csv, oldest first. An absent file reads as empty."""
    path = Path(blocks_dir) / INDEX_NAME
    if not path.is_file():
        return []
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def block_paths(blocks_dir: Path) -> List[Path]:
    """The block files in order. The zero-padded index makes the sort chronological."""
    blocks_dir = Path(blocks_dir)
    if not blocks_dir.is_dir():
        return []
    return sorted(blocks_dir.glob("block_*.log"))


def newest_block(blocks_dir: Path) -> Optional[Path]:
    blocks = block_paths(blocks_dir)
    return blocks[-1] if blocks else None


class LogSplitter:
    """Copy new bytes of an append-only log into bounded, indexed blocks.

    pump() is the whole interface. It is cheap when nothing has been written,
    it never raises for a reason the caller can act on, and it can be called at
    any interval: the block boundaries follow the byte count and the clock, not
    the call rate.
    """

    def __init__(self, source, blocks_dir=None, *,
                 max_bytes: Optional[int] = None,
                 max_seconds: Optional[float] = None):
        self.source = Path(source)
        self.blocks_dir = Path(blocks_dir) if blocks_dir else blocks_dir_for(self.source)
        self.max_bytes = block_bytes() if max_bytes is None else int(max_bytes)
        self.max_seconds = block_seconds() if max_seconds is None else float(max_seconds)

        self._offset = 0
        self._index = -1
        self._block: Optional[Path] = None
        self._opened_at: Optional[float] = None
        self._bytes = 0
        self._lines = 0
        self._first_eval: Optional[int] = None
        self._last_eval: Optional[int] = None
        self._carry = b""
        self._loaded = False

    # ------------------------------------------------------------------
    # State
    # ------------------------------------------------------------------
    @property
    def state_path(self) -> Path:
        return self.blocks_dir / STATE_NAME

    def _load_state(self) -> None:
        """Pick up where an earlier process stopped, if its state is usable."""
        if self._loaded:
            return
        self._loaded = True
        try:
            state = json.loads(self.state_path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            return
        try:
            self._offset = int(state["offset"])
            self._index = int(state["index"])
            name = state.get("block")
            self._block = self.blocks_dir / name if name else None
            self._opened_at = state.get("opened_at")
            self._bytes = int(state.get("bytes", 0))
            self._lines = int(state.get("lines", 0))
            self._first_eval = state.get("first_eval")
            self._last_eval = state.get("last_eval")
        except (KeyError, TypeError, ValueError):
            self._offset = 0
            self._index = -1
            self._block = None
        if self._block is not None and not self._block.is_file():
            # The block was removed under us. Start a new one rather than
            # recreating a file whose index row already says it closed.
            self._block = None
            self._opened_at = None
            self._bytes = 0
            self._lines = 0

    def _save_state(self) -> None:
        payload = {
            "source": str(self.source),
            "offset": self._offset,
            "index": self._index,
            "block": self._block.name if self._block else None,
            "opened_at": self._opened_at,
            "bytes": self._bytes,
            "lines": self._lines,
            "first_eval": self._first_eval,
            "last_eval": self._last_eval,
            "updated_at": _stamp(),
        }
        tmp = self.state_path.with_suffix(".json.tmp")
        try:
            tmp.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
            tmp.replace(self.state_path)
        except OSError:
            pass

    # ------------------------------------------------------------------
    # Blocks
    # ------------------------------------------------------------------
    def _open_block(self, note: str = "") -> None:
        self._index += 1
        self._opened_at = time.time()
        name = f"block_{self._index:04d}_{_stamp(self._opened_at)}.log"
        self._block = self.blocks_dir / name
        self._bytes = 0
        self._lines = 0
        self._first_eval = None
        self._last_eval = None
        header = (f"# block {self._index} opened "
                  f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(self._opened_at))}"
                  f" from {self.source.name} at byte {self._offset}")
        if note:
            header += f" ({note})"
        with self._block.open("wb") as handle:
            handle.write(header.encode("utf-8") + b"\n")
        try:
            (self.blocks_dir / LATEST_NAME).write_text(name + "\n", encoding="utf-8")
        except OSError:
            pass

    def _close_block(self, note: str = "") -> None:
        if self._block is None:
            return
        row = {
            "block": self._block.name,
            "opened_at": _stamp(self._opened_at),
            "closed_at": _stamp(),
            "bytes": self._bytes,
            "lines": self._lines,
            "first_eval": "" if self._first_eval is None else self._first_eval,
            "last_eval": "" if self._last_eval is None else self._last_eval,
            "note": note,
        }
        index_path = self.blocks_dir / INDEX_NAME
        try:
            new_file = not index_path.is_file()
            with index_path.open("a", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=INDEX_COLUMNS)
                if new_file:
                    writer.writeheader()
                writer.writerow(row)
        except OSError:
            pass
        self._block = None

    def _should_roll(self) -> bool:
        if self._block is None:
            return False
        if self._bytes >= self.max_bytes:
            return True
        return (self._opened_at is not None
                and (time.time() - self._opened_at) >= self.max_seconds)

    def _write_line(self, line: bytes) -> None:
        if self._block is None:
            self._open_block()
        try:
            with self._block.open("ab") as handle:
                handle.write(line)
        except OSError:
            return
        self._bytes += len(line)
        self._lines += 1
        found = EVAL_PATTERN.search(line)
        if found:
            eval_id = int(found.group(1))
            if self._first_eval is None:
                self._first_eval = eval_id
            self._last_eval = eval_id

    # ------------------------------------------------------------------
    # The pump
    # ------------------------------------------------------------------
    def pump(self) -> int:
        """Mirror whatever the source gained since the last call.

        Returns the number of bytes mirrored. Any filesystem error is reported
        once and swallowed: a mirror that fails must not stop a run.
        """
        try:
            return self._pump()
        except Exception as error:                       # noqa: BLE001
            print(f"[blocks] mirroring {self.source} failed: {error}",
                  file=sys.stderr, flush=True)
            return 0

    def _pump(self) -> int:
        if not self.source.is_file():
            return 0
        self.blocks_dir.mkdir(parents=True, exist_ok=True)
        self._load_state()

        size = self.source.stat().st_size
        if size < self._offset:
            # The log was replaced or truncated. Reading on from the old offset
            # would skip the whole new file.
            self._close_block("source truncated")
            self._offset = 0
            self._carry = b""
        if size == self._offset:
            if self._should_roll():
                self._close_block("time limit")
                self._save_state()
            return 0

        with self.source.open("rb") as handle:
            handle.seek(self._offset)
            chunk = handle.read(size - self._offset)
        self._offset += len(chunk)

        data = self._carry + chunk
        lines = data.splitlines(keepends=True)
        # A final line with no newline is still being written. Hold it back so
        # it is never split across two blocks.
        if lines and not lines[-1].endswith(b"\n"):
            self._carry = lines.pop()
        else:
            self._carry = b""

        for line in lines:
            if self._should_roll():
                self._close_block("size limit" if self._bytes >= self.max_bytes
                                  else "time limit")
            self._write_line(line)

        self._save_state()
        return len(chunk)

    def close(self) -> None:
        """Mirror the tail and close the open block, at the end of a run."""
        self.pump()
        if self._carry:
            self._write_line(self._carry + b"\n")
            self._carry = b""
        self._close_block("run ended")
        self._save_state()


# ----------------------------------------------------------------------
# Standalone use
# ----------------------------------------------------------------------
def _default_source(case: Optional[str]) -> Path:
    from run_config import resolve_results_dir
    base = Path(__file__).resolve().parents[1]
    return base / resolve_results_dir(case) / "logs" / "matlab_console.log"


def print_index(blocks_dir: Path) -> int:
    rows = read_index(blocks_dir)
    latest = newest_block(blocks_dir)
    if not rows and latest is None:
        print(f"[blocks] no blocks under {blocks_dir}.")
        return 1
    width = max([len(r["block"]) for r in rows] + [10])
    print(f"{'block'.ljust(width)}  {'opened':>15}  {'closed':>15}  "
          f"{'KiB':>8}  {'lines':>7}  evals")
    for row in rows:
        kib = int(row["bytes"] or 0) / 1024.0
        evals = f"{row['first_eval']}-{row['last_eval']}" if row["first_eval"] else ""
        print(f"{row['block'].ljust(width)}  {row['opened_at']:>15}  "
              f"{row['closed_at']:>15}  {kib:8.1f}  {row['lines']:>7}  {evals}")
    if latest is not None and latest.name not in {r["block"] for r in rows}:
        kib = latest.stat().st_size / 1024.0
        print(f"{latest.name.ljust(width)}  {'':>15}  {'(open)':>15}  {kib:8.1f}")
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("source", nargs="?", default=None,
                        help="the console log to mirror (default: the one under "
                             "the results folder that --case names)")
    parser.add_argument("--case", default=None,
                        help="the campaign: MFBO_RESULTS_DIR when set, else "
                             "results/<arm>_<case>_s<sobol_seed>")
    parser.add_argument("--follow", action="store_true",
                        help="keep mirroring until interrupted")
    parser.add_argument("--index", action="store_true",
                        help="print the block index and exit")
    parser.add_argument("--poll", type=float, default=2.0,
                        help="seconds between reads when following (default: 2.0)")
    args = parser.parse_args(sys.argv[1:] if argv is None else argv)

    source = Path(args.source) if args.source else _default_source(args.case)
    if args.index:
        return print_index(blocks_dir_for(source))

    splitter = LogSplitter(source)
    if not args.follow:
        moved = splitter.pump()
        print(f"[blocks] mirrored {moved} byte(s) into {splitter.blocks_dir}.")
        return 0

    print(f"[blocks] mirroring {source} into {splitter.blocks_dir}. Ctrl-C to stop.",
          flush=True)
    try:
        while True:
            splitter.pump()
            time.sleep(args.poll)
    except KeyboardInterrupt:
        splitter.close()
        print()
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
