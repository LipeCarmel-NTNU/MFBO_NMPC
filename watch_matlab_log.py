"""Print the MATLAB console log as it is written.

A headless MATLAB has no command window. run_supervised.py sends everything that
window would have shown to results/logs/matlab_console.log, and this script
prints that file line by line, so a second shell gives you the same view of the
server that the desktop gave you.

    python watch_matlab_log.py                 # last 40 lines, then follow
    python watch_matlab_log.py --lines 0       # follow only what arrives next
    python watch_matlab_log.py --all           # the whole file, then follow
    python watch_matlab_log.py --no-follow     # print and exit
    python watch_matlab_log.py results/logs/matlab_diary.log

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

# Same override as pipeline.matlab_interface.RESULTS_DIR: MFBO_RESULTS_DIR, default results.
DEFAULT_LOG = (Path(__file__).resolve().parent
               / os.environ.get("MFBO_RESULTS_DIR", "results") / "logs" / "matlab_console.log")


def parse(argv):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("path", nargs="?", default=str(DEFAULT_LOG),
                   help="the log to follow (default: results/logs/matlab_console.log)")
    p.add_argument("--lines", type=int, default=40,
                   help="lines of existing output to print before following (default: 40)")
    p.add_argument("--all", action="store_true",
                   help="print the whole file before following")
    p.add_argument("--no-follow", action="store_true",
                   help="print what is there and exit")
    p.add_argument("--poll", type=float, default=0.25,
                   help="seconds between reads (default: 0.25)")
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


def main(argv=None) -> int:
    args = parse(sys.argv[1:] if argv is None else argv)
    try:
        return follow(Path(args.path), args)
    except KeyboardInterrupt:
        print()
        return 130
    except BrokenPipeError:
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
