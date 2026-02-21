"""Load, document, and plot boxplot data exported by ``analyze_test_run_metrics.m``.

Figure description (focused subset used in Python workflows):
- Layout: 1x2 panels.
- Panel ``a``: settling-time boxplots by state.
  - x-axis states: ``V (L)``, ``X (g/L)``, ``S (g/L)``.
  - y-axis: settling time in hours.
  - includes a horizontal reference line at ``t = 8.09603 h``.
- Panel ``b`` (separate figure): prediction horizon ``N_p`` boxplot by case.

Boxplots follow a “scattered boxplot” style:
- Standard box (median, IQR, whiskers at 1.5 IQR)
- Overlaid raw data points with horizontal jitter
- Semi-transparent markers

Color policy:
- Use Wong (2011) Nature Methods palette entries via ``NATURE_HEX``.
- Use black for neutral scatter overlays.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman", "CMU Serif", "DejaVu Serif"],
        "mathtext.fontset": "cm",
        "axes.unicode_minus": False,
        "font.size": 12,
    }
)


PLOT_CONVENTION = {
    "panel_a": "Settling time by state (V, X, S)",
    "panel_b": "Prediction horizon N_p by case",
}

# Wong (2011) Nature Methods palette (hex forms for matplotlib use).
NATURE_HEX = {
    "Blue": "#0072B2",
    "BluishGreen": "#009E73",
    "ReddishPurple": "#CC79A7",
    "Vermillion": "#D55E00",
    "Orange": "#E69F00",
    "SkyBlue": "#56B4E9",
    "Yellow": "#F0E442",
    "Black": "#000000",
}

PLOT_CONVENTION.update(
    {
        "case_1_color_hex": NATURE_HEX["Blue"],
        "case_2_color_hex": NATURE_HEX["BluishGreen"],
        "guideline_color_hex": NATURE_HEX["ReddishPurple"],
        "scatter_color_hex": NATURE_HEX["Black"],
    }
)

REFERENCE_TIME_H = 8.09603
FINAL_PARETO_TIMESTAMPS = {
    "20260131_151035",
    "20260201_111557",
    "20260201_111928",
    "20260201_192807",
    "20260201_223337",
    "20260201_232106",
    "20260210_151703",
    "20260210_171107",
    "20260210_180826",
    "20260211_122653",
    "20260211_134235",
}


def load_boxplot_data(analysis_dir: str | Path | None = None) -> dict[str, pd.DataFrame]:
    """Load settling-time and final-error boxplot CSV exports."""
    base = Path(analysis_dir) if analysis_dir is not None else Path(__file__).resolve().parent
    project_root = base.parent
    settling_path = _first_existing(
        [base / "boxplot_settling_time_data.csv", project_root / "Result analysis" / "boxplot_settling_time_data.csv"]
    )
    final_err_path = _first_existing(
        [base / "boxplot_final_error_data.csv", project_root / "Result analysis" / "boxplot_final_error_data.csv"]
    )
    if settling_path is None or final_err_path is None:
        raise FileNotFoundError("Boxplot CSV files were not found. Run analyze_test_run_metrics.m first.")

    settling_df = pd.read_csv(settling_path)
    final_error_df = pd.read_csv(final_err_path)

    return {
        "settling_time": settling_df,
        "final_relative_error": final_error_df,
    }


def parse_np_from_summary(summary_txt_path: str | Path) -> pd.DataFrame:
    """Parse N_p rows from MATLAB text summary when CSV is unavailable."""
    path = Path(summary_txt_path)
    if not path.exists():
        return pd.DataFrame(columns=["run_label", "timestamp", "N_p"])

    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    start_idx = None
    for i, line in enumerate(lines):
        if "Controllers table (both case 1 and case 2, no removals):" in line:
            start_idx = i
            break
    if start_idx is None:
        return pd.DataFrame(columns=["run_label", "timestamp", "N_p"])

    records: list[dict[str, object]] = []
    for line in lines[start_idx + 1 :]:
        s = line.strip()
        if not s:
            continue
        if s.startswith("===") or s.startswith("NaN settling times") or s.startswith("Median settling times"):
            break
        if '"' not in s:
            continue

        parts = s.split('"')
        if len(parts) < 5:
            continue
        run_label = parts[1].strip()
        timestamp = parts[3].strip()
        tail = parts[4].strip()
        nums = tail.split()
        if len(nums) < 1:
            continue
        try:
            n_p = float(nums[0])
        except ValueError:
            continue
        records.append({"run_label": run_label, "timestamp": timestamp, "N_p": n_p})

    return pd.DataFrame.from_records(records)


def load_settling_and_np_data(analysis_dir: str | Path | None = None) -> dict[str, pd.DataFrame]:
    """Load settling-time data and N_p-by-case data with fallback logic."""
    data = load_boxplot_data(analysis_dir)
    base = Path(analysis_dir) if analysis_dir is not None else Path(__file__).resolve().parent
    project_root = base.parent

    summary_txt = _first_existing(
        [
            base / "analyze_test_run_metrics_summary.txt",
            project_root / "results" / "numerical results" / "analyze_test_run_metrics_summary.txt",
        ]
    )
    np_csv = _first_existing(
        [
            base / "boxplot_np_data.csv",
            project_root / "Result analysis" / "boxplot_np_data.csv",
        ]
    )

    if np_csv is not None:
        np_df = pd.read_csv(np_csv)
    elif summary_txt is not None:
        np_df = parse_np_from_summary(summary_txt)
    else:
        np_df = pd.DataFrame()

    if np_df.empty:
        np_df = build_np_from_results_csv(project_root)
    nc_df = build_nc_from_results_csv(project_root)
    return {"settling_time": data["settling_time"], "np_by_case": np_df, "nc_by_case": nc_df}


def build_np_from_results_csv(project_root: Path) -> pd.DataFrame:
    """Fallback: reconstruct N_p from run1/run2 results.csv and final Pareto timestamps."""
    records: list[dict[str, object]] = []
    for run_name, label in [("run1", "Case 1"), ("run2", "Case 2")]:
        csv_path = project_root / "results" / run_name / "results.csv"
        if not csv_path.exists():
            continue
        df = pd.read_csv(csv_path)
        if not {"timestamp", "theta_2", "theta_3"}.issubset(df.columns):
            continue

        df["timestamp"] = df["timestamp"].astype(str)
        df = df[df["timestamp"].isin(FINAL_PARETO_TIMESTAMPS)].copy()
        if df.empty:
            continue

        theta_p = np.rint(pd.to_numeric(df["theta_2"], errors="coerce"))
        theta_m = np.rint(pd.to_numeric(df["theta_3"], errors="coerce"))
        m = theta_m + 1
        n_p = theta_p + m

        for ts, val in zip(df["timestamp"].tolist(), n_p.tolist()):
            if pd.notna(val):
                records.append({"run_label": label, "timestamp": ts, "N_p": float(val)})

    return pd.DataFrame.from_records(records, columns=["run_label", "timestamp", "N_p"])


def build_nc_from_results_csv(project_root: Path) -> pd.DataFrame:
    """Reconstruct N_c from run1/run2 results.csv and final Pareto timestamps."""
    records: list[dict[str, object]] = []
    for run_name, label in [("run1", "Case 1"), ("run2", "Case 2")]:
        csv_path = project_root / "results" / run_name / "results.csv"
        if not csv_path.exists():
            continue
        df = pd.read_csv(csv_path)
        if not {"timestamp", "theta_3"}.issubset(df.columns):
            continue

        df["timestamp"] = df["timestamp"].astype(str)
        df = df[df["timestamp"].isin(FINAL_PARETO_TIMESTAMPS)].copy()
        if df.empty:
            continue

        theta_m = np.rint(pd.to_numeric(df["theta_3"], errors="coerce"))
        n_c = theta_m + 1

        for ts, val in zip(df["timestamp"].tolist(), n_c.tolist()):
            if pd.notna(val):
                records.append({"run_label": label, "timestamp": ts, "N_c": float(val)})

    return pd.DataFrame.from_records(records, columns=["run_label", "timestamp", "N_c"])


def _first_existing(paths: list[Path]) -> Path | None:
    """Return the first path that exists, else None."""
    for p in paths:
        if p.exists():
            return p
    return None


def scattered_boxplot(
    ax,
    data_groups: list[np.ndarray],
    labels: list[str],
    box_color: str,
    scatter_color: str,
    ylabel: str,
    show_scatter: bool = True,
):
    """Draw boxplot with optional jittered raw-point overlay."""
    positions = np.arange(1, len(data_groups) + 1)

    bp = ax.boxplot(
        data_groups,
        positions=positions,
        widths=0.5,
        patch_artist=True,
        showfliers=False,
        whis=1.5,
    )

    for patch in bp["boxes"]:
        patch.set(facecolor=box_color, alpha=0.25, edgecolor=box_color, linewidth=1, zorder=3)

    for median in bp["medians"]:
        median.set(color=box_color, linewidth=1.5, zorder=4)

    for whisker in bp["whiskers"]:
        whisker.set(color=box_color, linewidth=1, zorder=3)

    for cap in bp["caps"]:
        cap.set(color=box_color, linewidth=1, zorder=3)

    if show_scatter:
        # Overlay scattered raw data
        rng = np.random.default_rng(0)
        for i, y in enumerate(data_groups):
            x = np.full_like(y, positions[i], dtype=float)
            jitter = rng.uniform(-0.15, 0.15, size=len(y))
            ax.scatter(
                x + jitter,
                y,
                s=20,
                alpha=0.6,
                color=scatter_color,
                edgecolors="none",
                zorder=2,
            )

    ax.set_xticks(positions)
    ax.set_xticklabels(labels)
    ax.set_ylabel(ylabel)
    ax.grid(False)


def draw_guideline_behind_boxes(
    ax,
    y: float,
    positions: np.ndarray,
    box_width: float = 0.5,
    pad: float = 0.02,
    blocked_positions: list[float] | None = None,
    **line_kwargs,
):
    """Draw a horizontal guideline behind/around box bodies, not through them."""
    x_min, x_max = ax.get_xlim()
    half_w = box_width / 2.0 + pad
    blocked_src = positions if blocked_positions is None else blocked_positions
    blocked = sorted((float(p - half_w), float(p + half_w)) for p in blocked_src)

    cursor = x_min
    for left, right in blocked:
        if left > cursor:
            ax.hlines(y, cursor, left, **line_kwargs)
        cursor = max(cursor, right)
    if cursor < x_max:
        ax.hlines(y, cursor, x_max, **line_kwargs)


def plot_settling_time(df: pd.DataFrame):
    """Create settling-time boxplot by state with reference guideline."""
    fig, ax = plt.subplots(figsize=(6, 4))

    state_labels = ["V (L)", "X (g/L)", "S (g/L)"]
    data_groups: list[np.ndarray] = []

    # Expected CSV schema from MATLAB export is long-format:
    # run_label, timestamp, case_id, metric, state, state_label, value
    if {"state_label", "value"}.issubset(df.columns):
        for label in state_labels:
            vals = df.loc[df["state_label"] == label, "value"].dropna().to_numpy()
            data_groups.append(vals)
    elif {"state", "value"}.issubset(df.columns):
        state_ids = ["x1", "x2", "x3"]
        for sid in state_ids:
            vals = df.loc[df["state"] == sid, "value"].dropna().to_numpy()
            data_groups.append(vals)
    else:
        raise ValueError(
            "Unsupported settling dataframe schema. "
            "Expected long-format columns with state/state_label and value."
        )

    scattered_boxplot(
        ax=ax,
        data_groups=data_groups,
        labels=state_labels,
        box_color=PLOT_CONVENTION["case_1_color_hex"],
        scatter_color=PLOT_CONVENTION["scatter_color_hex"],
        ylabel="Settling time (h)",
        show_scatter=True,
    )

    draw_guideline_behind_boxes(
        ax=ax,
        y=REFERENCE_TIME_H,
        positions=np.arange(1, len(state_labels) + 1),
        blocked_positions=[3],
        box_width=0.5,
        color=PLOT_CONVENTION["guideline_color_hex"],
        linewidth=2.6,
        alpha=0.45,
        zorder=1,
    )

    fig.tight_layout()
    return fig


def plot_np_by_case(df: pd.DataFrame):
    """Create N_p-by-case boxplot with jittered raw-point overlay."""
    fig, ax = plt.subplots(figsize=(5, 4))

    if df is None or df.empty or "N_p" not in df.columns:
        ax.text(0.5, 0.5, "No $N_p$ data available", ha="center", va="center")
        ax.set_axis_off()
        fig.tight_layout()
        return fig

    df_local = df.copy()
    if "run_label" in df_local.columns:
        # Normalize common label variants from text parsing.
        df_local["run_label"] = (
            df_local["run_label"]
            .astype(str)
            .str.strip()
            .replace({"run1": "Case 1", "run2": "Case 2"})
        )
    else:
        df_local["run_label"] = "Case 1"

    df_local["N_p"] = pd.to_numeric(df_local["N_p"], errors="coerce")

    groups = []
    labels = []
    colors = []

    for case_label, color_key in [
        ("Case 1", "case_1_color_hex"),
        ("Case 2", "case_2_color_hex"),
    ]:
        subset = df_local[df_local["run_label"] == case_label]["N_p"].dropna().to_numpy()
        if len(subset) > 0:
            groups.append(subset)
            labels.append(case_label)
            colors.append(PLOT_CONVENTION[color_key])

    if len(groups) == 0:
        ax.text(0.5, 0.5, "No $N_p$ data available", ha="center", va="center")
        ax.set_axis_off()
        fig.tight_layout()
        return fig

    positions = np.arange(1, len(groups) + 1)

    bp = ax.boxplot(
        groups,
        positions=positions,
        widths=0.5,
        patch_artist=True,
        showfliers=False,
        whis=1.5,
    )

    for patch, color in zip(bp["boxes"], colors):
        patch.set(facecolor=color, alpha=0.25, edgecolor=color, linewidth=1)

    for median, color in zip(bp["medians"], colors):
        median.set(color=color, linewidth=1.5)

    rng = np.random.default_rng(0)
    for i, (y, color, case_label) in enumerate(zip(groups, colors, labels)):
        x = np.full_like(y, positions[i], dtype=float)
        jitter = rng.uniform(-0.15, 0.15, size=len(y))
        marker = "o"
        size = 36
        ax.scatter(
            x + jitter,
            y,
            s=size,
            marker=marker,
            alpha=0.6,
            color=PLOT_CONVENTION["scatter_color_hex"],
            edgecolors="none",
            zorder=2,
        )

    ax.set_xticks(positions)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Prediction horizon $N_p$")
    ax.set_ylim(bottom=0)
    ax.grid(False)

    fig.tight_layout()
    return fig


def plot_nc_by_case(df: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(5, 4))

    if df is None or df.empty or "N_c" not in df.columns:
        ax.text(0.5, 0.5, "No $N_c$ data available", ha="center", va="center")
        ax.set_axis_off()
        fig.tight_layout()
        return fig

    df_local = df.copy()
    if "run_label" in df_local.columns:
        df_local["run_label"] = (
            df_local["run_label"]
            .astype(str)
            .str.strip()
            .replace({"run1": "Case 1", "run2": "Case 2"})
        )
    else:
        df_local["run_label"] = "Case 1"

    df_local["N_c"] = pd.to_numeric(df_local["N_c"], errors="coerce")

    groups = []
    labels = []
    colors = []

    for case_label, color_key in [
        ("Case 1", "case_1_color_hex"),
        ("Case 2", "case_2_color_hex"),
    ]:
        subset = df_local[df_local["run_label"] == case_label]["N_c"].dropna().to_numpy()
        if len(subset) > 0:
            groups.append(subset)
            labels.append(case_label)
            colors.append(PLOT_CONVENTION[color_key])

    if len(groups) == 0:
        ax.text(0.5, 0.5, "No $N_c$ data available", ha="center", va="center")
        ax.set_axis_off()
        fig.tight_layout()
        return fig

    positions = np.arange(1, len(groups) + 1)

    bp = ax.boxplot(
        groups,
        positions=positions,
        widths=0.5,
        patch_artist=True,
        showfliers=False,
        whis=1.5,
    )

    for patch, color in zip(bp["boxes"], colors):
        patch.set(facecolor=color, alpha=0.25, edgecolor=color, linewidth=1)

    for median, color in zip(bp["medians"], colors):
        median.set(color=color, linewidth=1.5)

    rng = np.random.default_rng(0)
    for i, (y, color, case_label) in enumerate(zip(groups, colors, labels)):
        x = np.full_like(y, positions[i], dtype=float)
        jitter = rng.uniform(-0.15, 0.15, size=len(y))
        ax.scatter(
            x + jitter,
            y,
            s=36,
            marker="o",
            alpha=0.6,
            color=PLOT_CONVENTION["scatter_color_hex"],
            edgecolors="none",
            zorder=2,
        )

    ax.set_xticks(positions)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Control horizon $N_c$")
    ax.set_ylim(bottom=0)
    ax.grid(False)

    fig.tight_layout()
    return fig


def print_np_points_by_case(df: pd.DataFrame):
    """Print N_p values used in the boxplot, grouped by case label."""
    print("\nN_p points used in boxplot:")
    if df is None or df.empty or "N_p" not in df.columns:
        print("No N_p data available")
        return

    df_local = df.copy()
    if "run_label" in df_local.columns:
        df_local["run_label"] = (
            df_local["run_label"]
            .astype(str)
            .str.strip()
            .replace({"run1": "Case 1", "run2": "Case 2"})
        )
    else:
        df_local["run_label"] = "Case 1"

    df_local["N_p"] = pd.to_numeric(df_local["N_p"], errors="coerce")

    for case_label in ["Case 1", "Case 2"]:
        vals = df_local.loc[df_local["run_label"] == case_label, "N_p"].dropna().to_list()
        if vals:
            print(f"{case_label} ({len(vals)}): {vals}")
        else:
            print(f"{case_label} (0): []")


def print_nc_points_by_case(df: pd.DataFrame):
    print("\nN_c points used in boxplot:")
    if df is None or df.empty or "N_c" not in df.columns:
        print("No N_c data available")
        return

    df_local = df.copy()
    if "run_label" in df_local.columns:
        df_local["run_label"] = (
            df_local["run_label"]
            .astype(str)
            .str.strip()
            .replace({"run1": "Case 1", "run2": "Case 2"})
        )
    else:
        df_local["run_label"] = "Case 1"

    df_local["N_c"] = pd.to_numeric(df_local["N_c"], errors="coerce")

    for case_label in ["Case 1", "Case 2"]:
        vals = df_local.loc[df_local["run_label"] == case_label, "N_c"].dropna().to_list()
        if vals:
            print(f"{case_label} ({len(vals)}): {vals}")
        else:
            print(f"{case_label} (0): []")


if __name__ == "__main__":
    data = load_settling_and_np_data()

    settling_df = data["settling_time"]
    np_df = data["np_by_case"]
    nc_df = data["nc_by_case"]

    fig_a = plot_settling_time(settling_df)
    fig_b = plot_np_by_case(np_df)
    fig_c = plot_nc_by_case(nc_df)

    out_dir = Path(__file__).resolve().parent.parent / "results" / "graphical_results"
    out_dir.mkdir(parents=True, exist_ok=True)
    export_dpi = 600
    fig_a.savefig(out_dir / "python_settling_time_boxplot.png", dpi=300, bbox_inches="tight")
    fig_b.savefig(out_dir / "python_np_by_case_boxplot.png", dpi=300, bbox_inches="tight")
    fig_c.savefig(out_dir / "python_nc_by_case_boxplot.png", dpi=300, bbox_inches="tight")
    fig_a.savefig(out_dir / "python_settling_time_boxplot.pdf", dpi=export_dpi, bbox_inches="tight")
    fig_b.savefig(out_dir / "python_np_by_case_boxplot.pdf", dpi=export_dpi, bbox_inches="tight")
    fig_c.savefig(out_dir / "python_nc_by_case_boxplot.pdf", dpi=export_dpi, bbox_inches="tight")

    print_np_points_by_case(np_df)
    print_nc_points_by_case(nc_df)
    plt.show()
    print('Done')
