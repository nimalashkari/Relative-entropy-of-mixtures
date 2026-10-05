# Numerical report

The report TeX sources were moved from `report/` into this directory.

- Main document: `main.tex`
- Chapter sources: `sections/`
- Table sources: `tables/`
- Figures: `../figs/`
- Bibliography: included directly in `main.tex`; no BibTeX step required.

Compile from this directory with pdfLaTeX, repeating until cross-references settle:

```sh
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

When sharing or uploading the sources, include both `numerical/` and its sibling
`figs/`, preserving their relative layout.

The numerical data, generation scripts, and earlier build outputs remain in the
original `report/` directory. They are not required for compiling this document.
The generation scripts retain their original output paths; rerunning them does
not automatically update the relocated sources or figures.
