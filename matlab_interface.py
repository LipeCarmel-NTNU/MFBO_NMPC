"""Utilities for communicating theta with the MATLAB scripts."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable, Sequence
import time


BASE_DIR = Path(__file__).resolve().parent
LOCK_FILE = BASE_DIR / "matlab.lock"
THETA_FILE = BASE_DIR / "inbox" / "theta.txt"
RESULTS_FILE = BASE_DIR / "results" / "results.csv"


def read_results() -> tuple[list[str], list[float], list[float], list[float], list[float], list[list[float]]]:
    """Read optimization results ensuring lock coordination."""

    while LOCK_FILE.exists():
        time.sleep(5)

    if not RESULTS_FILE.exists():
        raise FileNotFoundError(f"Results file not found: {RESULTS_FILE}")

    timestamp: list[str] = []
    sse: list[float] = []
    ssd_u: list[float] = []
    cost: list[float] = []
    runtime: list[float] = []
    theta_matrix: list[list[float]] = []

    with RESULTS_FILE.open("r", encoding="ascii", newline="") as fh:
        reader = csv.DictReader(fh)
        required_columns = [
            "timestamp",
            "SSE",
            "SSdU",
            "J",
            "runtime_s",
            *[f"theta_{i}" for i in range(1, 13)],
        ]
        missing = [column for column in required_columns if column not in reader.fieldnames]
        if missing:
            raise ValueError(f"Results CSV missing columns: {', '.join(missing)}")

        for row in reader:
            timestamp.append(row["timestamp"])
            sse.append(float(row["SSE"]))
            ssd_u.append(float(row["SSdU"]))
            cost.append(float(row["J"]))
            runtime.append(float(row["runtime_s"]))
            theta_row = [float(row[f"theta_{i}"]) for i in range(1, 13)]
            theta_matrix.append(theta_row)

    return timestamp, sse, ssd_u, cost, runtime, theta_matrix

def send_theta(theta: Sequence[float] | Iterable[float]) -> None:
    """Validate theta and forward it to MATLAB via ``write_theta``."""

    values = _validate_theta(theta)
    write_theta(values)


def _is_integer_value(value: object) -> bool:
    """Return True when value behaves like a non-boolean integer."""

    if isinstance(value, bool):
        return False
    if isinstance(value, int):
        return True
    return isinstance(value, float) and value.is_integer()


def _is_positive_real(value: object) -> bool:
    """Return True when value is a numeric type greater than zero."""

    if isinstance(value, bool):
        return False
    if isinstance(value, (int, float)):
        return value > 0
    return False


def _validate_theta(theta: Sequence[float] | Iterable[float]) -> list[float]:
    """Validate and normalize theta parameters."""

    values = list(theta)
    expected_length = 12
    if len(values) != expected_length:
        raise ValueError(
            f"theta must contain {expected_length} elements, received {len(values)}"
        )

    # Validate f: must be in [0, 1]
    f = values[0]
    if not _is_positive_real(f) or not (0 <= f <= 1):
        raise ValueError(
            f"theta[0]: f must be a real number in [0, 1], got {f}"
        )
    values[0] = float(f)

    # Validate theta_p: non-negative integer
    theta_p = values[1]
    if not _is_integer_value(theta_p) or int(theta_p) < 0:
        raise ValueError(
            f"theta[1]: theta_p must be a non-negative integer, got {theta_p}"
        )
    values[1] = int(round(float(theta_p)))

    # Validate theta_m: non-negative integer
    theta_m = values[2]
    if not _is_integer_value(theta_m) or int(theta_m) < 0:
        raise ValueError(
            f"theta[2]: theta_m must be a non-negative integer, got {theta_m}"
        )
    values[2] = int(round(float(theta_m)))

    # Validate q values: positive reals
    for i in range(3):
        q_value = values[3 + i]
        if not _is_positive_real(q_value):
            raise ValueError(
                f"theta[{3 + i}]: q[{i}] must be a positive real number, got {q_value}"
            )

    # Validate r^(u) values: positive reals
    for i in range(3):
        r_u_value = values[6 + i]
        if not _is_positive_real(r_u_value):
            raise ValueError(
                f"theta[{6 + i}]: r^(u)[{i}] must be a positive real number, got {r_u_value}"
            )

    # Validate r^(Delta u) values: positive reals
    for i in range(3):
        r_du_value = values[9 + i]
        if not _is_positive_real(r_du_value):
            raise ValueError(
                f"theta[{9 + i}]: r^(Delta u)[{i}] must be a positive real number, got {r_du_value}"
            )

    return values



def write_theta(theta: Sequence[float]) -> None:
    """Write theta values to ``inbox/theta.txt`` after lock negotiation."""

    while LOCK_FILE.exists():
        time.sleep(5)

    THETA_FILE.parent.mkdir(parents=True, exist_ok=True)
    with THETA_FILE.open("w", encoding="ascii") as fh:
        fh.write(" ".join(str(value) for value in theta))
        fh.write("\n")


