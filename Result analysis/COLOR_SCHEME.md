# Result analysis — figure color scheme

Reference for every plotting script in `Result analysis/`. The conventions below
were extracted from `old/main_pareto.m` (the manuscript figure hub) and the
helpers in `dependencies/plot_utils/`. New scripts (e.g. `surrogate_analysis.m`)
must follow this scheme.

## 1. Categorical palette — Wong (2011) Nature Methods

Source: `dependencies/plot_utils/nature_methods_colors.m`
(Wong, B. *Points of view: Color blindness.* Nat Methods 8, 441 (2011),
doi:10.1038/nmeth.1618). Colorblind-safe.

`plotColors = nature_methods_colors(3)` returns, in this order:

| Row | Name          | RGB (0–255)     | Hex       | Role in figures                          |
|-----|---------------|-----------------|-----------|------------------------------------------|
| 1   | Blue          | (0, 114, 178)   | `#0072B2` | **Case 1** (lines, markers `o`)          |
| 2   | BluishGreen   | (0, 158, 115)   | `#009E73` | **Case 2** (lines, markers `^`); also "promoted" points |
| 3   | ReddishPurple | (204, 121, 167) | `#CC79A7` | Accent: Pareto frontier, guide curves, fidelity-z overlays |

> Note: two stale comments in `old/main_pareto.m` (lines ~60 and ~78) claim the
> triple is "Blue, Vermillion, Orange". That is wrong — the function's actual
> order is Blue, BluishGreen, ReddishPurple (the comment near line 897 is the
> correct one). Trust the function, not those comments.

Remaining Wong colors, available via `nature_methods_colors()` (struct) when
more series are needed: Vermillion `#D55E00`, Orange `#E69F00`,
SkyBlue `#56B4E9`, Yellow `#F0E442`.

## 2. Sequential colormap — navia (Scientific Colour Maps, Crameri)

Source: `dependencies/plot_utils/navia/navia.mat` (licence and README in the
same folder). Loaded via a local `load_navia_colormap(n)` helper that
interpolates the stored map to `n` levels and clips to [0, 1].

**Truncated and reversed, and there are now two of them.** `ra_context` calls
`ra_seq_colormap(name, ctx.seqCrFloor)` with `seqCrFloor = 3.0`:

| Field | Map | Used by |
|---|---|---|
| `ctx.seqMapZ` | navia | fidelity `z` on the objective-space markers |
| `ctx.seqMapVintage` | batlow | vintage index on the surrogate phi(z) curves |

`ctx.seqMap` remains as an alias of `seqMapZ` for generators that still read it.

A perceptually uniform sequential map is uniform in *lightness*, which is what
drives one end of it to paper-white: navia ends at `#FCF4D9`, 1.10:1 against
white. That is right for a filled field, which is what these maps are drawn for,
and useless for discrete markers and thin lines on white paper. `ra_seq_colormap`
therefore orients the map dark-first, cuts the pale end where contrast against
white falls through `seqCrFloor`, and reverses it so the HIGH value takes the
dark end. Cut this way every map ends near L\* 61, so the choice between them is
a choice of hue path, not of how much was trimmed. Lowering `seqCrFloor`
lengthens the ramp and weakens its faintest colour.

Two maps rather than one because `z` and the vintage index are different
quantities in different figures; a shared ramp would read as one scale carried
across the figure set.

Maps are stored in `dependencies/plot_utils/colormaps/seq_maps.mat`. The README
beside it carries provenance and the citation each map needs — **batlow needs the
same Crameri citation as navia**, so adding it changes nothing in the
bibliography. `Result analysis/pick_seq_colormap.m` renders the candidates side
by side.

**Colour axis.** `caxis` spans the pooled data, never a fixed floor:
`gen_data_frontier` stores `zLo`/`zHi` as the min and max of `z` over every case
in the run, so the panels of one figure share a scale and equal colour means
equal `z` across them. The old `min(0.5, min(z))` floor spent two thirds of the
ramp on an empty range and flattened every marker to one colour.

Used for any *continuous* quantity mapped to color:

- fidelity `z` in the SSE-vs-SSdU scatter (`caxis [0.5, 1]`, colorbar labeled
  `$z$ (dimensionless)`)
- any ordered progression (e.g. surrogate vintage index) — one navia sample
  per level, colorbar labeled accordingly

## 3. Auxiliary colors

| Color                   | Value                     | Role                                            |
|-------------------------|---------------------------|--------------------------------------------------|
| Black                   | `k`                       | Raw optimization-sample scatter (filled), DOE cutoff `xline` (dashed, LineWidth 2), axis colors |
| Benchmark orange        | (242, 133, 34) `#F28522`  | Benchmark controller markers; matches `dependencies/plot_utils/good_colors.m` `C.orange` |
| Reference cloud blue    | `[0.60 0.82 0.98]`        | Re-evaluated reference sample cloud (light, background) |
| White                   | `w`                       | Figure background, open-marker faces on frontier continuum |

`good_colors.m` (red `#FF1F5B`, green `#00CD6C`, blue `#009ADE`, purple
`#AF58BA`, yellow `#FFC61E`, orange `#F28522`) is a secondary palette; in the
results analysis only its orange is used, for benchmarks.

## 4. Styling conventions (non-color, but part of the look)

- `figure("Color", "w")`; `grid off`, `box off` on every axes.
- LaTeX everywhere: `set(groot, "defaultTextInterpreter", "latex")` plus the
  axes-tick and legend equivalents.
- Font size **20** on axes, labels, and colorbars (`set_font_size` /
  `set(ax, "FontSize", fontSize)`).
- Panel labels: bold lowercase letters as titles, `$\mathbf{a}$`,
  `$\mathbf{b}$`, ..., with `TitleHorizontalAlignment = "left"`.
- Tick labels formatted with `format_tick(Xdec, Ydec)` (powers of 10 on log
  axes).
- Multi-panel: `tiledlayout(..., "Padding", "compact", "TileSpacing", "compact")`,
  side-by-side panels ~1200x460 px.
- Export: PNG at 300 dpi (`exportgraphics`) plus vector PDF with the same stem
  (`save_figure`), standardized by a `save_plot_outputs(fig, pngPath, fontSize, W, H)`
  helper. Figures land in `results/graphical_results/`.

## 5. Quick recipe for a new script

```matlab
addpath(genpath(fullfile(repo_root, "dependencies")));
plotColors = nature_methods_colors(3); % Blue, BluishGreen, ReddishPurple
seqMapZ       = ra_seq_colormap("navia",  3.0); % markers: truncated, high = dark
seqMapVintage = ra_seq_colormap("batlow", 3.0); % surrogate curves
% Case 1 -> plotColors(1,:), marker "o"
% Case 2 -> plotColors(2,:), marker "^"
% accent / frontier -> plotColors(3,:)
```
