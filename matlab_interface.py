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


def _is_in_range(value: object, min_val: float, max_val: float) -> bool:
    """Return True when value is a numeric type in [min_val, max_val]."""

    if isinstance(value, bool):
        return False
    if isinstance(value, (int, float)):
        return min_val <= value <= max_val
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
    if not _is_in_range(f, 0, 1):
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

    # Validate q values: in [-3, 3]
    for i in range(3):
        q_value = values[3 + i]
        if not _is_in_range(q_value, -3, 3):
            raise ValueError(
                f"theta[{3 + i}]: q[{i}] must be in [-3, 3], got {q_value}"
            )

    # Validate r^(u) values: intentionally disabled.
    # We allow very negative exponents (e.g., -1000) so 10^r_u is effectively zero in MATLAB.
    # for i in range(3):
    #     r_u_value = values[6 + i]
    #     if not _is_in_range(r_u_value, -3, 3):
    #         raise ValueError(
    #             f"theta[{6 + i}]: r^(u)[{i}] must be in [-3, 3], got {r_u_value}"
    #         )

    # Validate r^(Delta u) values: in [-3, 3]
    for i in range(3):
        r_du_value = values[9 + i]
        if not _is_in_range(r_du_value, -3, 3):
            raise ValueError(
                f"theta[{9 + i}]: r^(Delta u)[{i}] must be in [-3, 3], got {r_du_value}"
            )

    return values


def send_theta(theta: Sequence[float] | Iterable[float]) -> None:
    """Validate theta and forward it to MATLAB via ``write_theta``."""

    values = _validate_theta(theta)
    write_theta(values)


def write_theta(theta: Sequence[float]) -> None:
    """Write theta values to ``inbox/theta.txt`` after lock negotiation."""

    while LOCK_FILE.exists():
        time.sleep(5)

    THETA_FILE.parent.mkdir(parents=True, exist_ok=True)
    with THETA_FILE.open("w", encoding="ascii") as fh:
        fh.write(" ".join(str(value) for value in theta))
        fh.write("\n")


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


if __name__ == "__main__":
    # Tuning vector theta components
    f = 0.01        # Fraction/Fidelity parameter in (0, 1)
    theta_p = 17    # Prediction horizon offset, integer in [0, 60]
    theta_m = 2     # Control horizon offset, integer in [0, 30]

    # State cost weights (in log10 scale, maps to Q = diag(10^q))
    q = [1.3010299956639813, 0, 0]
    
    # Input cost weights (in log10 scale, maps to R_u = diag(10^r^(u)))
    r_u = [0.3010299956639812, 0.3010299956639812, 0]
    
    # Input rate cost weights (in log10 scale, maps to R_Delta_u = diag(10^r^(Delta u)))
    r_delta_u = [2, 2, 1]
    
    # Assemble theta vector
    dummy_theta = [f, theta_p, theta_m] + q + r_u + r_delta_u
    
    send_theta(dummy_theta)


