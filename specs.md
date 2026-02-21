# Figure Formatting Specs

This file defines the current project-wide plotting conventions used in analysis scripts (for example, `resultssandbox.m` and `surrogate.m`).

## Global style
- Use LaTeX interpreters for text, ticks, and legends.
- Default font size: `14`.
- Disable axis grid and box (`grid off`, `box off`) unless a specific figure requires otherwise.
- Use consistent figure dimensions via `set_fig_size(...)`.
- For panel labels used as titles (for example, `\textbf{a}` and `\textbf{b}`), left-align the title inside each axes.
- In vertically stacked panels with matching and explicit x-axes, omit the x-label on the top panel and keep it only on the bottom panel.

## Colors
- Use `good_colors(n)` from `dependencies/plot_utils`.
- For paired line plots, use:
  - color 1 for measured/data curve,
  - color 2 for fitted/model curve.

## Lines and markers
- Prefer line width `2.0` for primary curves.
- Use markers only when they add information (for example, sampled points or highlighted Pareto points).

## Tick labels
- Use `format_tick(X_decimals, Y_decimals)` for fixed decimal formatting.
- For log axes, `format_tick` uses scientific labels (`10^k` or `m x 10^k`).

## Export policy
- Save figures as PNG (raster) and PDF (vector) through shared helpers.
- Place exported graphics under `results/graphical_results/`.

## Pareto plots
- Use log-log axes when plotting `J_{TV}` vs `J_{track}`.
- Keep guideline/continuum curves visually behind markers.
- Use consistent axis limits when comparing cases across figures.
