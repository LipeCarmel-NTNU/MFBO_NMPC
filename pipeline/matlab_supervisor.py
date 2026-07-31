"""Own the MATLAB server process and relaunch it when it stops serving.

The driver blocks on a results row for up to eval_timeout_s. A MATLAB process
that has died produces no row, so without supervision an interruption costs the
whole timeout and then a manual restart of both halves. This module makes the
MATLAB process a child of the driver's process, which turns that situation into
a signal the driver can act on: the exit of a known pid.

Two conditions trigger a relaunch.

  process exit   The pid is gone. This needs no timer.
  no lock        A request is pending, the pid is alive, and matlab.lock has
                 been absent for wedge_timeout_s. serve_requests holds that
                 lock for as long as it works on a request, so its absence
                 while a request is unanswered means the work is not running.
                 This covers a process that stopped serving without exiting.

The relaunch count is unbounded. A crash that recurs is answered again, because
stopping a run that would otherwise have continued costs more than the cores
that a repeated relaunch spends. A launch that never reaches the serve loop is
paced by an increasing delay so that a fault at startup, an unreachable license
server for instance, does not fill the log with thousands of attempts.

MATLAB runs headless. -nodisplay removes the X11 connection, so no dropped ssh
session can take the process down, and the console output that the desktop would
show is written to results/logs/matlab_console.log instead. watch_matlab_log.py
follows that file.

One detail of the MATLAB side has to be compensated here. main_initialization.m
deletes inbox/theta.txt on startup so that a request from an earlier run is not
the first thing it serves. A relaunch during the design phase therefore removes
the request that the driver is waiting for. This module reads the request before
it kills the process and writes it back once the new process is serving.
"""

from __future__ import annotations

import csv
import os
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import List, Optional, Tuple

from pipeline import matlab_interface as mi

BASE_DIR = Path(__file__).resolve().parents[1]
LOG_DIR = BASE_DIR / "results" / "logs"
CONSOLE_LOG = LOG_DIR / "matlab_console.log"
RESTART_LOG = LOG_DIR / "restarts.csv"

# serve_requests prints this line once it is polling the inbox. It prints it
# after the startup deletions, which makes it the point from which a request is
# safe to publish. Both entry points reach it through the same function.
READY_MARKER = "Serving requests from"

ENTRY_FOR_PHASE = {mi.PHASE_DOE: "main_initialization", mi.PHASE_OPT: "main_BO"}

RESTART_COLUMNS = ["timestamp", "attempt", "reason", "exit_code", "entry",
                   "pending_eval_id", "republished"]


class MatlabSupervisor:
    """A MATLAB server process that is restarted for as long as it is needed.

    The instance is not thread-safe beyond its own use: start() and stop() run
    on the calling thread, and the monitor runs on one thread of its own.
    """

    def __init__(self, entry: str = "main_initialization", *,
                 matlab: str = "matlab",
                 ready_timeout_s: float = 900.0,
                 wedge_timeout_s: float = 600.0,
                 poll_s: float = 5.0,
                 kill_grace_s: float = 20.0,
                 backoff_max_s: float = 300.0,
                 diary: bool = False,
                 console_log: Path = CONSOLE_LOG,
                 restart_log: Path = RESTART_LOG):
        self.entry = entry
        self.matlab = matlab
        self.ready_timeout_s = ready_timeout_s
        self.wedge_timeout_s = wedge_timeout_s
        self.poll_s = poll_s
        self.kill_grace_s = kill_grace_s
        self.backoff_max_s = backoff_max_s
        self.diary = diary
        self.console_log = Path(console_log)
        self.restart_log = Path(restart_log)

        self._proc: Optional[subprocess.Popen] = None
        self._thread: Optional[threading.Thread] = None
        self._stop = threading.Event()
        self._launches = 0
        self._launch_offset = 0
        self._startup_failures = 0
        self._lock_absent_since: Optional[float] = None
        self.restarts = 0

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def start(self) -> None:
        """Launch MATLAB and return once it is serving requests.

        The caller must not publish a request before this returns. A request
        written during the startup of main_initialization would be deleted by
        it, and the driver would then wait out eval_timeout_s for a row that
        nobody is going to write.
        """
        while not self._stop.is_set():
            self._launch(self.entry)
            if self._await_ready():
                self._startup_failures = 0
                return
            # The launch did not reach the serve loop. Whatever stopped it is
            # not specific to a request, so the next attempt is delayed.
            self._startup_failures += 1
            delay = min(5.0 * 2 ** (self._startup_failures - 1), self.backoff_max_s)
            self._log(f"the launch did not reach the serve loop "
                      f"({self._startup_failures} in a row). Retrying in {delay:.0f} s.")
            self._tail_console(20)
            self._sleep(delay)

    def watch(self) -> None:
        """Start the monitor on its own thread."""
        if self._thread is not None:
            return
        self._thread = threading.Thread(target=self._monitor, name="matlab-supervisor",
                                        daemon=True)
        self._thread.start()

    def stop(self) -> None:
        """End the monitor and the MATLAB process."""
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=self.poll_s * 2 + 5.0)
            self._thread = None
        if self.is_alive():
            self._log("stopping MATLAB.")
            self._terminate()
        mi.LOCK_FILE.unlink(missing_ok=True)

    def is_alive(self) -> bool:
        return self._proc is not None and self._proc.poll() is None

    # ------------------------------------------------------------------
    # Monitor
    # ------------------------------------------------------------------

    def _monitor(self) -> None:
        while not self._stop.wait(self.poll_s):
            reason, exit_code = self._failure()
            if reason is None:
                continue
            if self._stop.is_set():
                return

            pending = self._pending_request()
            saved = self._read_theta_raw() if pending else None
            entry = ENTRY_FOR_PHASE.get(pending[1], self.entry) if pending else self.entry

            self.restarts += 1
            self._log(f"restart {self.restarts}: {reason}. "
                      f"Relaunching {entry}"
                      + (f" for pending evaluation {pending[0]}." if pending else "."))

            self._terminate()
            # acquire_lock warns about a lock it did not create, and the driver
            # treats one as a busy server. Neither is true after a kill.
            mi.LOCK_FILE.unlink(missing_ok=True)
            self._lock_absent_since = None

            republished = False
            while not self._stop.is_set():
                self._launch(entry)
                if self._await_ready():
                    self._startup_failures = 0
                    republished = self._republish(saved, pending)
                    break
                self._startup_failures += 1
                delay = min(5.0 * 2 ** (self._startup_failures - 1), self.backoff_max_s)
                self._log(f"the relaunch did not reach the serve loop "
                          f"({self._startup_failures} in a row). Retrying in {delay:.0f} s.")
                self._tail_console(20)
                self._sleep(delay)

            self._record_restart(reason, exit_code, entry, pending, republished)

    def _failure(self) -> Tuple[Optional[str], Optional[int]]:
        """Name the condition that calls for a relaunch, or return None."""
        if self._proc is None:
            return "no MATLAB process", None

        code = self._proc.poll()
        if code is not None:
            return f"the MATLAB process exited with code {code}", code

        # A pending request with no lock means the request is not being worked
        # on. The lock is absent for a moment between evaluations as well, so
        # the condition has to hold for wedge_timeout_s. The margin also covers
        # the handover from main_initialization to main_BO, which happens inside
        # this same process and leaves the lock absent while it runs.
        if self._pending_request() is None:
            self._lock_absent_since = None
            return None, None

        if mi.LOCK_FILE.exists():
            self._lock_absent_since = None
            return None, None

        now = time.time()
        if self._lock_absent_since is None:
            self._lock_absent_since = now
            return None, None

        absent_s = now - self._lock_absent_since
        if absent_s > self.wedge_timeout_s:
            return (f"the process is alive but {mi.LOCK_FILE.name} has been absent "
                    f"for {absent_s / 60:.1f} min with a request pending"), None
        return None, None

    # ------------------------------------------------------------------
    # The MATLAB process
    # ------------------------------------------------------------------

    def _launch(self, entry: str) -> None:
        self.console_log.parent.mkdir(parents=True, exist_ok=True)
        self._launches += 1

        statement = entry
        if self.diary:
            # diary writes each line as it is produced. It is available for the
            # case where the redirection below turns out to be buffered.
            diary_path = str(self.console_log.with_name("matlab_diary.log")).replace("'", "''")
            statement = f"diary('{diary_path}'); {entry}"

        cmd = [self.matlab] + self._headless_flags() + ["-batch", statement]

        banner = (f"\n{'=' * 70}\n"
                  f"launch {self._launches}: {' '.join(cmd)}\n"
                  f"{time.strftime('%Y-%m-%d %H:%M:%S')}\n"
                  f"{'=' * 70}\n")
        self._launch_offset = self._console_size()
        with self.console_log.open("ab") as fh:
            fh.write(banner.encode("utf-8"))
            fh.flush()

        # The log is handed to MATLAB as a file descriptor, so its output goes
        # to the file directly and no buffer of this process sits in between.
        log_fh = self.console_log.open("ab")
        try:
            self._proc = subprocess.Popen(
                cmd, cwd=str(BASE_DIR), stdout=log_fh, stderr=subprocess.STDOUT,
                stdin=subprocess.DEVNULL, **self._process_group_kwargs())
        except FileNotFoundError:
            self._proc = None
            self._log(f"'{self.matlab}' is not on PATH. On a module system, load the "
                      f"MATLAB module before starting this command.")
        finally:
            log_fh.close()

        if self._proc is not None:
            self._log(f"MATLAB launched as pid {self._proc.pid} "
                      f"({entry}, console {self.console_log}).")

    def _headless_flags(self) -> List[str]:
        if os.name == "nt":
            return ["-nosplash", "-nodesktop"]
        return ["-nodisplay", "-nodesktop", "-nosplash"]

    def _process_group_kwargs(self) -> dict:
        """Put MATLAB in a group of its own so the pool workers can be signalled."""
        if os.name == "nt":
            return {"creationflags": subprocess.CREATE_NEW_PROCESS_GROUP}
        return {"start_new_session": True}

    def _await_ready(self) -> bool:
        """Wait for the serve loop, or return False if the launch failed."""
        if self._proc is None:
            return False
        deadline = time.time() + self.ready_timeout_s
        start_size = self._console_size()
        while not self._stop.is_set():
            if self._console_contains(READY_MARKER, since=start_size):
                self._log("MATLAB is serving requests.")
                return True
            if self._proc.poll() is not None:
                self._log(f"MATLAB exited with code {self._proc.returncode} "
                          f"before it reached the serve loop.")
                return False
            if time.time() > deadline:
                self._log(f"MATLAB did not reach the serve loop within "
                          f"{self.ready_timeout_s:.0f} s.")
                self._terminate()
                return False
            time.sleep(2.0)
        return False

    def _terminate(self) -> None:
        if self._proc is None or self._proc.poll() is not None:
            self._proc = None if self._proc is None else self._proc
            return
        try:
            if os.name == "nt":
                subprocess.run(["taskkill", "/F", "/T", "/PID", str(self._proc.pid)],
                               capture_output=True, check=False)
            else:
                # The pool workers are children of the MATLAB client and share
                # its process group, so the group is what has to be signalled.
                os.killpg(os.getpgid(self._proc.pid), signal.SIGTERM)
            self._proc.wait(timeout=self.kill_grace_s)
        except (ProcessLookupError, PermissionError, OSError):
            pass
        except subprocess.TimeoutExpired:
            try:
                if os.name != "nt":
                    os.killpg(os.getpgid(self._proc.pid), signal.SIGKILL)
                self._proc.wait(timeout=10.0)
            except (ProcessLookupError, PermissionError, OSError,
                    subprocess.TimeoutExpired):
                pass

    # ------------------------------------------------------------------
    # The pending request
    # ------------------------------------------------------------------

    def _read_theta_raw(self) -> Optional[str]:
        try:
            return mi.THETA_FILE.read_text(encoding="ascii")
        except OSError:
            return None

    def _pending_request(self) -> Optional[Tuple[int, int]]:
        """Return (eval_id, phase_code) of a request with no row yet.

        A request that already has a results row or a failure row is answered,
        and the driver has moved on or is about to.
        """
        raw = self._read_theta_raw()
        if raw is None:
            return None
        fields = raw.split()
        if len(fields) != 2 + mi.THETA_LEN:
            return None
        try:
            eval_id = int(float(fields[0]))
            phase_code = int(float(fields[1]))
        except ValueError:
            return None
        if phase_code not in ENTRY_FOR_PHASE:
            return None

        phase = "init" if phase_code == mi.PHASE_DOE else "bo"
        try:
            if any(r["eval_id"] == eval_id for r in mi.read_results(mi.results_file(phase))):
                return None
            if any(f["eval_id"] == eval_id
                   for f in mi.read_failures(mi.failures_file(phase))):
                return None
        except (OSError, ValueError):
            return None
        return eval_id, phase_code

    def _republish(self, saved: Optional[str],
                   pending: Optional[Tuple[int, int]]) -> bool:
        """Write back a request that the restarted server deleted on startup.

        The write is skipped when the file is present, which means either that
        the server kept it or that the driver has since published a newer one,
        and when the saved request has been answered in the meantime.
        """
        if saved is None or pending is None:
            return False
        if mi.THETA_FILE.exists():
            return False
        if self._pending_request() is not None:
            return False

        eval_id = pending[0]
        phase = "init" if pending[1] == mi.PHASE_DOE else "bo"
        try:
            if any(r["eval_id"] == eval_id for r in mi.read_results(mi.results_file(phase))):
                return False
        except (OSError, ValueError):
            pass

        tmp = mi.THETA_FILE.with_name("theta.supervisor.tmp")
        try:
            mi.THETA_FILE.parent.mkdir(parents=True, exist_ok=True)
            with tmp.open("w", encoding="ascii") as fh:
                fh.write(saved)
                fh.flush()
                os.fsync(fh.fileno())
            tmp.replace(mi.THETA_FILE)
        except OSError as exc:
            self._log(f"could not write the pending request back to "
                      f"{mi.THETA_FILE}: {exc}")
            return False
        self._log(f"the relaunched server deleted the inbox on startup, so "
                  f"evaluation {eval_id} was published again.")
        return True

    # ------------------------------------------------------------------
    # Records and output
    # ------------------------------------------------------------------

    def _record_restart(self, reason: str, exit_code: Optional[int], entry: str,
                        pending: Optional[Tuple[int, int]], republished: bool) -> None:
        self.restart_log.parent.mkdir(parents=True, exist_ok=True)
        new = not self.restart_log.exists()
        try:
            with self.restart_log.open("a", encoding="ascii", newline="") as fh:
                writer = csv.writer(fh)
                if new:
                    writer.writerow(RESTART_COLUMNS)
                writer.writerow([
                    time.strftime("%Y%m%d_%H%M%S"), self.restarts, reason,
                    "" if exit_code is None else exit_code, entry,
                    "" if pending is None else pending[0], int(republished),
                ])
        except OSError as exc:
            self._log(f"could not append to {self.restart_log}: {exc}")

    def _console_size(self) -> int:
        try:
            return self.console_log.stat().st_size
        except OSError:
            return 0

    def _console_contains(self, marker: str, since: int = 0) -> bool:
        """Look for a marker in the part of the console log written since an offset.

        The offset matters because the log is appended across launches. Without
        it, the marker printed by an earlier launch would satisfy the barrier of
        the current one immediately.
        """
        try:
            with self.console_log.open("rb") as fh:
                fh.seek(since)
                return marker.encode("utf-8") in fh.read()
        except OSError:
            return False

    def _tail_console(self, n_lines: int) -> None:
        """Print the end of what the current launch wrote.

        The offset keeps the output of earlier launches out of it. Without it,
        every failed attempt would reprint the whole history of the run.
        """
        try:
            with self.console_log.open("rb") as fh:
                fh.seek(self._launch_offset)
                lines = fh.read().splitlines()[-n_lines:]
        except OSError:
            return
        if not lines:
            return
        print(f"[matlab] last {len(lines)} line(s) of {self.console_log}:", flush=True)
        for line in lines:
            print("  " + line.decode("utf-8", "replace"), flush=True)

    def _sleep(self, seconds: float) -> None:
        self._stop.wait(seconds)

    def _log(self, message: str) -> None:
        print(f"[matlab] {message}", flush=True)
        sys.stdout.flush()
