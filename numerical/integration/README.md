# Numerical integration into the theory manuscript

The main entry point is `../../Relative Entropy Mixture.tex`.
Compile from `manuscript_theory/`, including this directory and `figs/`.

All numerical prose, tables, captions and appendix methods are now written
directly in that main TeX file; it does not input the historical TeX fragments
in this directory. Edit the main source for manuscript changes. Figure labels
follow its notation: relative entropy S, complement fraction epsilon=1-x,
length l, ordered states (Psi,Phi), and fit models M_1 and M_2. CSV state and
model keys remain internal to the renderer and are mapped to these labels.

- Section 2.2: direct primary-state comparison and a twelve-direction fit table.
- Section 4.1: compact-boson small- and large-interval results, two fit tables, four figures.
- Section 4.2: Ising small- and large-interval results, two fit tables, four figures.
- Appendix L: lattice states, full sampling grid, reduced matrices, transition/Gram entropy, stable entropy remainders, precision checks, and fitting methods.

The integration uses the adopted lengths 6--10 and windows (small Ising,
small boson, large Ising, large boson) = (.02, .1, .002, .005).
It does not use archived length extensions or replace these windows with
the historical .0047 large-boson selection.

The four condensed mixture tables copy values from the corresponding
`../tables/{ising,boson}_{small,large}_parameters.tex` files.
The primary table copies `../tables/{ising,boson}_primary_all_fits.tex`.
The eight manuscript figures are `../../figs/*_{small,large}_fit_window_{value,residual}.pdf`.
They show only fitted observations and curves within the adopted windows.
All lengths share one circular marker and one data color; residual colors
identify fitting models, not subsystem lengths. The grid is copied from
`../tables/extended_grid_definition.tex` with distinct labels. Source tables,
original full-range figure PDFs, physical observations, and fits are unchanged.

Regenerate these figure variants from the stored CSVs, from the project root:

```sh
MPLCONFIGDIR=/tmp/relative_entropy_mpl /Library/Frameworks/Python.framework/Versions/3.14/bin/python3 manuscript_theory/numerical/integration/plot_fit_windows.py --audit artifacts/task_manuscript_plot_windows/plot_audit.json
```

The renderer reads `report/data/scaling_observations.csv`,
`scaling_fit_summary.csv`, and `scaling_predictions.csv`, and checks fitted
coefficients against `f2_summary_fits.csv`. It does not recompute fits.

The numerical methods are condensed from `../sections/01_construction.tex`,
`01_ising_fermions.tex`, `05_small_x_extension.tex`, `11_complement.tex`,
`14_larger_method.tex`, `14_scaling_fits.tex`, and `15_window_selection.tex`.
The precision diagnostics come from `../tables/scaling_numerical_checks.tex`;
window counts and sensitivity statements come from `common_window_counts.tex`
and `density_conclusions.tex`. Bibliographic records reproduce the supplied
numerical report's method references, using distinct citation keys.

For a PDF without editorial margin-label stamps, run from `manuscript_theory/`:

```sh
pdflatex -interaction=nonstopmode -halt-on-error -jobname=relative_entropy_mixture '\PassOptionsToPackage{final}{showlabels}\input{Relative Entropy Mixture.tex}'
bibtex relative_entropy_mixture
pdflatex -interaction=nonstopmode -halt-on-error -jobname=relative_entropy_mixture '\PassOptionsToPackage{final}{showlabels}\input{Relative Entropy Mixture.tex}'
pdflatex -interaction=nonstopmode -halt-on-error -jobname=relative_entropy_mixture '\PassOptionsToPackage{final}{showlabels}\input{Relative Entropy Mixture.tex}'
```

The source retains its original `showlabels` setting. Omitting the
`PassOptionsToPackage` command produces the author's annotated view.
The delivered `Relative Entropy Mixture.pdf` uses the clean view.

The original draft has two unresolved analytical references (`fnexpression`
and `eq:F_n_1_minus_x_mixed`), pre-existing overfull equations, and undefined
BibTeX month strings. They are outside this numerical integration. All
numerical labels and citations resolve. The integration diff, original
source snapshot, and archived superseded numerical discussion are retained
in `artifacts/task_manuscript_numerics/` at the project root.
