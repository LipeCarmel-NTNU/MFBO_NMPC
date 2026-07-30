"""Copy the console output of a run to a file.

The console is the only place where the progress of a run appears in order:
the design points, the proposals, the refits, the wall times and any warning
that a library raises. A terminal keeps a limited scrollback, and a run lasts
hours, so that record has to reach a file.

This module writes both streams. Text still goes to the terminal, and a copy
goes to results/logs/<phase>_<timestamp>.log. The copy is line-buffered and
flushed, so the file is complete up to the last line even if the process is
killed.

A new file is created for each start. An interrupted run therefore leaves its
log in place, and the restart writes a second file next to it.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import List, TextIO


class _Tee:
    """Write to several streams at once.

    Only the methods that print and the warnings module use are needed. The
    file stream is flushed on every write, because the value of this file is
    highest when the process dies without closing it.
    """

    def __init__(self, streams: List[TextIO], flush_always: List[bool]):
        self._streams = streams
        self._flush_always = flush_always

    def write(self, text: str) -> int:
        for stream, always in zip(self._streams, self._flush_always):
            try:
                stream.write(text)
                if always:
                    stream.flush()
            except (ValueError, OSError):
                pass          # a closed or full stream must not stop the run
        return len(text)

    def flush(self) -> None:
        for stream in self._streams:
            try:
                stream.flush()
            except (ValueError, OSError):
                pass

    def isatty(self) -> bool:
        return getattr(self._streams[0], "isatty", lambda: False)()

    @property
    def encoding(self) -> str:
        return getattr(self._streams[0], "encoding", "utf-8")


class ConsoleLog:
    """Context manager that mirrors stdout and stderr into a file.

    Use it around a whole phase:

        with ConsoleLog(root, "bo") as log:
            ...
            print(f"the record is in {log.path}")
    """

    def __init__(self, root: Path, phase: str, subdir: str = "logs"):
        self.dir = Path(root) / subdir
        self.phase = phase
        stamp = time.strftime("%Y%m%d_%H%M%S")
        self.path = self.dir / f"{phase}_{stamp}.log"
        self._file: TextIO | None = None
        self._stdout = None
        self._stderr = None

    def __enter__(self) -> "ConsoleLog":
        self.dir.mkdir(parents=True, exist_ok=True)
        self._file = self.path.open("w", encoding="utf-8", buffering=1)
        self._file.write(f"# {self.phase} phase, started "
                         f"{time.strftime('%Y-%m-%dT%H:%M:%S%z')}\n")
        self._file.flush()

        self._stdout, self._stderr = sys.stdout, sys.stderr
        sys.stdout = _Tee([self._stdout, self._file], [False, True])
        sys.stderr = _Tee([self._stderr, self._file], [False, True])
        print(f"[log] console output is copied to {self.path}")
        return self

    def __exit__(self, exc_type, exc, tb) -> bool:
        if exc_type is not None:
            # The traceback that follows goes to the file as well, because
            # sys.stderr is still the tee at this point.
            print(f"[log] the phase ended on {exc_type.__name__}")
        if self._stdout is not None:
            sys.stdout, sys.stderr = self._stdout, self._stderr
        if self._file is not None:
            self._file.write(f"# ended {time.strftime('%Y-%m-%dT%H:%M:%S%z')}\n")
            self._file.flush()
            self._file.close()
            self._file = None
        return False          # never swallow an exception
