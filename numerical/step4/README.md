# Complete step-4 equal-mixture results

The active manuscript is `Relative Entropy Mixture.tex`. Its four equal-mixture numerical tables and eight figures use `fits.csv`, copied exactly from the verified full-grid analysis. The common fitting windows are `0 < x <= 0.02` and `0 < epsilon <= 0.005` for both CFTs, with lengths 6 through 10 and every permitted system size divisible by four.

| Regime | Model | Observations per direction | Distinct sizes |
|---|---|---:|---:|
| Small | Ising | 9744 | 1974 |
| Small | Compact boson | 19984 | 4022 |
| Large | Ising | 8245 | 1749 |
| Large | Compact boson | 18485 | 3797 |

The figure generator is `numerical/integration/plot_uniform_step4.py`. Run it from the surrounding research workspace with its installed NumPy/Matplotlib environment. It reads `cluster_uniform_step4/results/full_observations.csv`, checks the pinned data hash, and renders every selected observation. Scatter layers are rasterized at 350 dpi inside the PDF; axes, text, and fitted curves remain vector graphics. No observations are thinned or averaged.

The source dataset contains exactly 365148 verified measurements, including 181980 large-interval deficits. SHA256: `32c3a3d54eee116a7d5e17d48d6074fa92c41d6bcf7e59dbb6eb889412ead5e9`. `figure_audit.json` records the selected counts, sizes, fractions, coordinate hashes, and residual ranges for all 48 panels. The manuscript can compile from the supplied static figures without the source dataset or Python environment.

The black CFT curves contain the leading term and any complete known next term: small-Ising vacuum/energy sixth-order terms and large-boson sixth- and third-order terms. Other channels show the leading term only. Fitted large-interval coefficients describe `log(2)-D`; their signs are therefore opposite to those in the entropy expansion. Relative errors are signed percentages computed before rounding.

The general-y figures retain their separate 14-size data and stated windows. Historical scripts and unused numerical fragments are retained for provenance; they do not supply the active equal-mixture tables or figures. The manuscript's updated text explicitly distinguishes fit agreement from a continuum extrapolation and does not claim a qualifying common large window.
