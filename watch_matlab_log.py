"""Print the MATLAB console log as it is written.

A headless MATLAB has no command window. run_supervised.py sends everything that
window would have shown to <results folder>/logs/matlab_console.log, and this
script prints that file line by line, so a second shell gives you the same view
of the server that the desktop gave you. --case names the campaign, and the
folder resolves as the driver resolves it: MFBO_RESULTS_DIR when it is set,
else results/<arm>_<case>_s<sobol_seed>.

    python watch_matlab_log.py --case baseline              # last 40 lines, then follow
    python watch_matlab_log.py --case baseline --lines 0    # follow only what arrives next
    python watch_matlab_log.py --case baseline --all        # the whole file, then follow
    python watch_matlab_log.py --case baseline --no-follow  # print and exit
    python watch_matlab_log.py results/SF_baseline_s123/logs/matlab_diary.log

The raw log is one appended file for the whole campaign, so reading it from the
top is slow by the end of a run. run_supervised.py also mirrors it into bounded
blocks under logs/console_blocks/ (see pipeline/log_blocks.py). Two options read
that view instead:

    python watch_matlab_log.py --case baseline --list    # the block index, newest last
    python watch_matlab_log.py --case baseline --blocks  # follow the newest block

--blocks follows the block being written and moves to the next one when it
opens, so the file open at any moment is at most one block long. --list prints
which evaluation ids each closed block holds, which is how you find the block
that covers the evaluation you care about.

Every line is flushed as it is printed, so the output is complete up to the last
line even when this script is piped into another command or killed.

The script waits for the file when it does not exist yet, and it reopens the file
when it is replaced or truncated. It therefore survives a restart of the run and
does not need to be started in any particular order.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

from run_config import CASES, resolve_results_dir

BASE_DIR = Path(__file__).resolve().parent


def parse(argv):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("path", nargs="?", default=None,
                   help="the log to follow (default: <results folder>/logs/"
                        "matlab_console.log, the folder that --case names)")
    p.add_argument("--case", default=None, choices=sorted(CASES),
                   help="the campaign whose log to follow, resolved as the driver "
                        "resolves it: MFBO_RESULTS_DIR when set, else "
                        "results/<arm>_<case>_s<sobol_seed>")
    p.add_argument("--lines", type=int, default=40,
                   help="lines of existing output to print before following (default: 40)")
    p.add_argument("--all", action="store_true",
                   help="print the whole file before following")
    p.add_argument("--no-follow", action="store_true",
                   help="print what is there and exit")
    p.add_argument("--poll", type=float, default=0.25,
                   help="seconds between reads (default: 0.25)")
    p.add_argument("--blocks", action="store_true",
                   help="follow the newest console block instead of the raw log")
    p.add_argument("--list", action="store_true", dest="list_blocks",
                   help="print the console block index and exit")
    return p.parse_args(argv)


def emit(chunk: bytes) -> None:
    """Write bytes through to stdout and flush.

    The log is written by MATLAB, so it can hold a byte sequence that is not
    valid UTF-8, and one bad byte must not end the follow. Decoding with
    "replace" keeps every other character.
    """
    sys.stdout.write(chunk.decode("utf-8", "replace"))
    sys.stdout.flush()


def start_offset(handle, args) -> int:
    """Where to begin reading, given how much history was asked for."""
    size = handle.seek(0, os.SEEK_END)
    if args.all:
        return 0
    if args.lines <= 0:
        return size

    # Read backwards in blocks until the requested number of line breaks is in
    # hand. A log of any size is handled without reading all of it.
    block = 8192
    offset = size
    newlines = 0
    while offset > 0 and newlines <= args.lines:
        step = min(block, offset)
        offset -= step
        handle.seek(offset)
        newlines += handle.read(step).count(b"\n")
    handle.seek(offset)
    lines = handle.read(size - offset).splitlines(keepends=True)
    if offset > 0 and lines:
        # The scan stopped inside a line, so the first entry is a fragment.
        lines = lines[1:]
    keep = lines[-args.lines:]
    return size - sum(len(line) for line in keep)


def follow(path: Path, args) -> int:
    announced = False
    while not path.exists():
        if args.no_follow:
            print(f"[watch] {path} does not exist.", file=sys.stderr)
            return 1
        if not announced:
            print(f"[watch] waiting for {path} to appear ...", flush=True)
            announced = True
        time.sleep(1.0)

    handle = path.open("rb")
    try:
        stat = os.fstat(handle.fileno())
        position = start_offset(handle, args)
        handle.seek(position)

        while True:
            chunk = handle.read(65536)
            if chunk:
                emit(chunk)
                continue

            if args.no_follow:
                return 0

            # A new file at the same path, or a file that shrank, means the log
            # was replaced or truncated. Reading on from the old offset would
            # either miss the new content or return nothing at all.
            try:
                current = path.stat()
            except OSError:
                time.sleep(args.poll)
                continue

            replaced = (current.st_ino, current.st_dev) != (stat.st_ino, stat.st_dev)
            truncated = current.st_size < handle.tell()
            if replaced or truncated:
                print(f"\n[watch] {path} was "
                      f"{'replaced' if replaced else 'truncated'}. Reopening.",
                      flush=True)
                handle.close()
                handle = path.open("rb")
                stat = os.fstat(handle.fileno())
                continue

            time.sleep(args.poll)
    finally:
        handle.close()


def blocks_dir_for(log_path: Path) -> Path:
    """Where the blocks of a given log live. Kept in step with pipeline.log_blocks."""
    return log_path.parent / "console_blocks"


def sorted_blocks(blocks_dir: Path):
    """The block files in order. The zero-padded index makes the sort chronological."""
    if not blocks_dir.is_dir():
        return []
    return sorted(blocks_dir.glob("block_*.log"))


def print_block_index(blocks_dir: Path) -> int:
    """Print index.csv, plus the block still being written."""
    index = blocks_dir / "index.csv"
    blocks = sorted_blocks(blocks_dir)
    if not index.is_file() and not blocks:
        print(f"[watch] no console blocks under {blocks_dir}.", file=sys.stderr)
        print("[watch] they are written by run_supervised.py, or by "
              "python -m pipeline.log_blocks --follow.", file=sys.stderr)
        return 1
    closed = set()
    if index.is_file():
        import csv
        with index.open("r", newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        width = max([len(r["block"]) for r in rows] + [10])
        print(f"{'block'.ljust(width)}  {'opened':>15}  {'KiB':>8}  {'lines':>7}  evals")
        for row in rows:
            closed.add(row["block"])
            kib = int(row["bytes"] or 0) / 1024.0
            span = (f"{row['first_eval']}-{row['last_eval']}"
                    if row["first_eval"] else "")
            print(f"{row['block'].ljust(width)}  {row['opened_at']:>15}  "
                  f"{kib:8.1f}  {row['lines']:>7}  {span}")
    for block in blocks:
        if block.name not in closed:
            print(f"{block.name}  (open, {block.stat().st_size / 1024.0:.1f} KiB)")
    return 0


def follow_blocks(blocks_dir: Path, args) -> int:
    """Follow the newest block, moving on when the next one opens."""
    announced = False
    while not sorted_blocks(blocks_dir):
        if args.no_follow:
            print(f"[watch] no console blocks under {blocks_dir}.", file=sys.stderr)
            return 1
        if not announced:
            print(f"[watch] waiting for a block under {blocks_dir} ...", flush=True)
            announced = True
        time.sleep(1.0)

    current = sorted_blocks(blocks_dir)[-1]
    print(f"[watch] following {current.name}", flush=True)
    handle = current.open("rb")
    try:
        handle.seek(start_offset(handle, args))
        while True:
            chunk = handle.read(65536)
            if chunk:
                emit(chunk)
                continue
            if args.no_follow:
                return 0

            newest = sorted_blocks(blocks_dir)[-1]
            if newest != current:
                # Drain the block that just closed before moving on, so the last
                # lines written to it are not lost.
                emit(handle.read())
                handle.close()
                current = newest
                print(f"\n[watch] block rolled. Following {current.name}", flush=True)
                handle = current.open("rb")
                continue
            time.sleep(args.poll)
    finally:
        handle.close()


def main(argv=None) -> int:
    args = parse(sys.argv[1:] if argv is None else argv)
    path = (Path(args.path) if args.path else
            BASE_DIR / resolve_results_dir(args.case) / "logs" / "matlab_console.log")
    try:
        if args.list_blocks:
            return print_block_index(blocks_dir_for(path))
        if args.blocks:
            return follow_blocks(blocks_dir_for(path), args)
        return follow(path, args)
    except KeyboardInterrupt:
        print()
        return 130
    except BrokenPipeError:
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
