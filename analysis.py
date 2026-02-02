from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from matlab_interface import read_results


DEFAULT_PLOT_PATH = Path("results") / "sse_vs_ssdu.png"
RUNTIME_PLOT_PATH = Path("results") / "runtime_vs_iteration.png"


THETA_COLUMN_RENAME = {
    "theta_1": "f",
    "theta_2": "theta_p",
    "theta_3": "theta_m",
    "theta_4": "q1_log10",
    "theta_5": "q2_log10",
    "theta_6": "q3_log10",
    "theta_7": "r_u1_log10",
    "theta_8": "r_u2_log10",
    "theta_9": "r_u3_log10",
    "theta_10": "r_du1_log10",
    "theta_11": "r_du2_log10",
    "theta_12": "r_du3_log10",
}


SERIES_LABELS = ["x1", "x2", "x3"]


def _decode_theta(theta_matrix: Iterable[Iterable[float]]) -> pd.DataFrame:
    """Return a DataFrame with named theta components and derived weights."""

    theta_df = pd.DataFrame(
        theta_matrix,
        columns=[f"theta_{i}" for i in range(1, 13)],
        dtype=float,
    ).rename(columns=THETA_COLUMN_RENAME)

    # Convert integral parameters back to ints after validation in MATLAB layer.
    theta_df["theta_p"] = theta_df["theta_p"].round().astype(int)
    theta_df["theta_m"] = theta_df["theta_m"].round().astype(int)

    # Translate log-scaled weights (base-10) to actual diagonal entries.
    for idx, label in enumerate(SERIES_LABELS, start=1):
        theta_df[f"Q_{label}"] = 10.0 ** theta_df[f"q{idx}_log10"]
        theta_df[f"R_u_{label}"] = 10.0 ** theta_df[f"r_u{idx}_log10"]
        theta_df[f"R_du_{label}"] = 10.0 ** theta_df[f"r_du{idx}_log10"]

    return theta_df


def build_results_dataframe() -> pd.DataFrame:
    """Load the optimization results CSV and attach decoded theta columns."""

    timestamp, sse, ssd_u, cost, runtime, theta = read_results()
    base_df = pd.DataFrame(
        {
            "timestamp": pd.to_datetime(timestamp, format="%Y%m%d_%H%M%S"),
            "SSE": sse,
            "SSdU": ssd_u,
            "J": cost,
            "runtime_s": runtime,
        }
    )
    theta_df = _decode_theta(theta)
    df = pd.concat([base_df, theta_df], axis=1)
    df["iteration"] = range(1, len(df) + 1)
    df["runtime_min"] = df["runtime_s"] / 60.0
    df["p"] = df["theta_p"] + df["theta_m"]
    df["m"] = df["theta_m"]
    return df


def plot_sse_vs_ssdu(df: pd.DataFrame):
    """Create a scatter plot comparing SSE and SSdU."""

    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(9, 6))

    scatter = sns.scatterplot(
        data=df,
        x="SSE",
        y="SSdU",
        hue="J",
        size="theta_p",
        palette="viridis",
        sizes=(40, 300),
        ax=ax,
        linewidth=0.7,
        edgecolor="black",
        alpha=0.9,
    )

    ax.set_title("SSE vs SSdU")
    ax.set_xlabel("SSE (State Tracking Error)")
    ax.set_ylabel("SSdU (Control Effort)")
    ax.grid(True, which="major", linestyle="--", linewidth=0.5)
    sns.despine(ax=ax)
    scatter.legend(title="Cost / Horizon", loc="best", frameon=True)
    fig.tight_layout()

    return fig, ax


def plot_runtime_vs_iteration(df: pd.DataFrame):
    """Plot runtime (minutes) versus iteration index."""

    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(9, 5.5))

    ax.plot(
        df["iteration"],
        df["runtime_min"],
        marker="o",
        linewidth=2.0,
        markersize=6,
        color="#2a7f62",
    )

    ax.axvline(
        20,
        color="#d7263d",
        linestyle="--",
        linewidth=1.3,
        label="Optimization start (iter 20)",
    )

    ax.set_xlabel("Iteration")
    ax.set_ylabel("Runtime (minutes)")
    ax.set_title("Iteration Runtime")
    ax.grid(True, which="major", linestyle="--", linewidth=0.5)
    sns.despine(ax=ax)
    ax.legend(loc="best", frameon=True)
    fig.tight_layout()

    return fig, ax


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Analyze MFBO results")
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_PLOT_PATH,
        help="Path to save the SSE vs SSdU scatter plot",
    )
    parser.add_argument(
        "--runtime-output",
        type=Path,
        default=RUNTIME_PLOT_PATH,
        help="Path to save the runtime vs iteration plot",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display the plot interactively after saving",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    df = build_results_dataframe()
    sse_fig, _ = plot_sse_vs_ssdu(df)
    runtime_fig, _ = plot_runtime_vs_iteration(df)

    output_path = args.output
    runtime_output = args.runtime_output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    runtime_output.parent.mkdir(parents=True, exist_ok=True)
    sse_fig.savefig(output_path, dpi=300, bbox_inches="tight")
    runtime_fig.savefig(runtime_output, dpi=300, bbox_inches="tight")

    if args.show:
        plt.show()
    else:
        plt.close(sse_fig)
        plt.close(runtime_fig)

    print(f"Saved SSE vs SSdU scatter plot to {output_path}")
    print(f"Saved runtime vs iteration plot to {runtime_output}")


if __name__ == "__main__":
    main()
