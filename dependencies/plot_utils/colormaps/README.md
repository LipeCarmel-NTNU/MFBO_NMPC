# seq_maps.mat — sequential colormaps for the chooser

One 256x3 sRGB matrix per variable, sampled uniformly over each map's range.
Written once so `Result analysis/pick_seq_colormap.m` can compare candidates
without a network call. `navia` here is the same map already vendored under
`dependencies/plot_utils/navia/`, resampled the same way.

## Provenance and licence

- `navia batlow lajolla davos oslo bamako imola nuuk lapaz turku acton bilbao
  grayC hawaii buda devon tokyo` — Scientific Colour Maps, Fabio Crameri.
  Licensed CC BY 4.0. Cite: Crameri, F. (2018), Scientific colour maps,
  Zenodo, doi:10.5281/zenodo.1243862; and Crameri, F., G. E. Shephard and
  P. J. Heron (2020), The misuse of colour in science communication,
  Nature Communications 11, 5444, doi:10.1038/s41467-020-19160-7.
  Sampled via the `cmcrameri` Python package, which redistributes Crameri's
  published map data unchanged.
- `cividis` — Nunez, J. R., C. R. Anderton and R. S. Renslow (2018), PLoS ONE
  13(7), e0199239. Public domain, as shipped with Matplotlib.
- `viridis` — Smith, N. J. and S. van der Walt. CC0, as shipped with Matplotlib.
- `YlGnBu`, `Blues` — ColorBrewer, Cynthia Brewer, Mark Harrower and
  The Pennsylvania State University. Apache 2.0, as shipped with Matplotlib.

Whatever ends up in the paper figures should be cited in the caption or the
methods; the Crameri reference above is the one the current figures need.
