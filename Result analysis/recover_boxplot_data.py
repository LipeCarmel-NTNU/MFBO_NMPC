"""Load and document boxplot data exported by ``analyze_test_run_metrics.m``.

Figure description (focused subset used in Python workflows):
- Layout: 1x2 panels.
- Panel ``a``: settling-time boxplots by state.
  - x-axis states: ``V (L)``, ``X (g/L)``, ``S (g/L)``.
  - y-axis: settling time in hours.
  - includes a horizontal reference line at ``t = 8.09603 h``.
- Panel ``b`` (separate figure): prediction horizon ``N_p`` boxplot by case.

Color convention used across result-analysis figures:
- Case 1 / panel ``a``: color index 1 from the active palette
  (hex ``#FF1F5B``, RGB ``(255, 31, 91)``).
- Case 2 / panel ``b``: color index 2 from the active palette
  (hex ``#009ADE``, RGB ``(0, 154, 222)``).
- Guideline / fidelity overlays use color index 3
  (hex ``#AF58BA``, RGB ``(175, 88, 186)``).

Note:
- In ``resultssandbox.m``, colors are generated via ``good_colors(4)`` and
  then the 4th row is removed. This leaves ``[red, blue, purple]``.
"""

from __future__ import annotations

from pathlib import Path
import pandas as pd

PLOT_CONVENTION = {
    "panel_a": "Settling time by state (V, X, S)",
    "panel_b": "Final relative error by state (V, X, S)",
    "case_1_color_hex": "#FF1F5B",
    "case_2_color_hex": "#009ADE",
    "guideline_color_hex": "#AF58BA",
}


def load_boxplot_data(analysis_dir: str | Path | None = None) -> dict[str, pd.DataFrame]:
    """Return exported boxplot data as pandas DataFrames.

    Parameters
    ----------
    analysis_dir:
        Folder containing:
        - boxplot_settling_time_data.csv
        - boxplot_final_error_data.csv
        Defaults to the directory of this script.
    """
    base = Path(analysis_dir) if analysis_dir is not None else Path(__file__).resolve().parent

    settling_path = base / "boxplot_settling_time_data.csv"
    final_err_path = base / "boxplot_final_error_data.csv"

    settling_df = pd.read_csv(settling_path)
    final_error_df = pd.read_csv(final_err_path)

    return {
        "settling_time": settling_df,
        "final_relative_error": final_error_df,
    }


def load_settling_and_np_data(analysis_dir: str | Path | None = None) -> dict[str, pd.DataFrame]:
    """Return only settling-time and N_p dataframes.

    - settling_time: from ``boxplot_settling_time_data.csv``
    - np_by_case: inferred from per-controller table in
      ``analyze_test_run_metrics_summary.txt`` (both-case, no removals)
    """
    data = load_boxplot_data(analysis_dir)
    base = Path(analysis_dir) if analysis_dir is not None else Path(__file__).resolve().parent
    summary_txt = base / "analyze_test_run_metrics_summary.txt"
    np_df = parse_np_from_summary(summary_txt)
    return {"settling_time": data["settling_time"], "np_by_case": np_df}


def parse_np_from_summary(summary_txt_path: str | Path) -> pd.DataFrame:
    """Parse N_p by case from analyze_test_run_metrics_summary.txt.

    Expected source section:
    'Controllers table (both case 1 and case 2, no removals):'
    """
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

        # Minimal robust parse: run_label and timestamp are quoted, then numeric columns.
        # Example:
        # "Case 1"  "2026..."  1  1  ...
        parts = s.split('"')
        if len(parts) < 5:
            continue
        run_label = parts[1].strip()
        timestamp = parts[3].strip()
        tail = parts[4].strip()
        nums = tail.split()
        if len(nums) < 2:
            continue
        try:
            n_p = float(nums[0])
        except ValueError:
            continue
        records.append({"run_label": run_label, "timestamp": timestamp, "N_p": n_p})

    return pd.DataFrame.from_records(records)


def get_plot_convention() -> dict[str, str]:
    """Return textual plotting convention used in the MATLAB analysis figures."""
    return dict(PLOT_CONVENTION)


if __name__ == "__main__":
    data = load_settling_and_np_data()
    print("Loaded DataFrames (focused set):")
    for name, df in data.items():
        print(f"- {name}: {df.shape[0]} rows x {df.shape[1]} cols")
    print("\nPlot convention:")
    for k, v in get_plot_convention().items():
        print(f"- {k}: {v}")
