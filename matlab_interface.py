"""Utilities for communicating theta with the MATLAB scripts."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence
import time


BASE_DIR = Path(__file__).resolve().parent
LOCK_FILE = BASE_DIR / "matlab.lock"
THETA_FILE = BASE_DIR / "inbox" / "theta.txt"


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


def send_theta(theta: Sequence[float] | Iterable[float]) -> None:
    """Validate theta and forward it to MATLAB via ``write_theta``."""

    values = list(theta)
    expected_length = 12
    if len(values) != expected_length:
        raise ValueError(
            f"theta must contain {expected_length} elements, received {len(values)}"
        )

    iter_max = values[0]
    if not _is_integer_value(iter_max) or not (1 <= int(iter_max) <= 300):
        raise ValueError(
            f"theta[0]: Iter_max must be an integer in [1, 300], got {iter_max}"
        )
    values[0] = int(round(float(iter_max)))

    theta_p = values[1]
    if not _is_integer_value(theta_p) or int(theta_p) < 0:
        raise ValueError(
            f"theta[1]: theta_p must be a non-negative integer, got {theta_p}"
        )
    values[1] = int(round(float(theta_p)))

    theta_m = values[2]
    if not _is_integer_value(theta_m) or int(theta_m) < 0:
        raise ValueError(
            f"theta[2]: theta_m must be a non-negative integer, got {theta_m}"
        )
    values[2] = int(round(float(theta_m)))

    for i in range(3):
        q_value = values[3 + i]
        if not _is_positive_real(q_value):
            raise ValueError(
                f"theta[{3 + i}]: q[{i}] must be a positive real number, got {q_value}"
            )

    for i in range(3):
        r_u_value = values[6 + i]
        if not _is_positive_real(r_u_value):
            raise ValueError(
                f"theta[{6 + i}]: r^(u)[{i}] must be a positive real number, got {r_u_value}"
            )

    for i in range(3):
        r_du_value = values[9 + i]
        if not _is_positive_real(r_du_value):
            raise ValueError(
                f"theta[{9 + i}]: r^(Delta u)[{i}] must be a positive real number, got {r_du_value}"
            )

    write_theta(values)


def write_theta(theta: Sequence[float]) -> None:
    """Write theta values to ``inbox/theta.txt`` after lock negotiation."""

    while LOCK_FILE.exists():
        time.sleep(5)

    THETA_FILE.parent.mkdir(parents=True, exist_ok=True)
    with THETA_FILE.open("w", encoding="ascii") as fh:
        fh.write(" ".join(str(value) for value in theta))
        fh.write("\n")
