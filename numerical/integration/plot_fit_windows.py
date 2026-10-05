"""Render the eight manuscript plots from stored observations and fits.

No fitting or density-matrix calculation is performed. Every plotted
observation is selected by the adopted fit window. Subsystem length is
not encoded in marker shape or color; residual colors identify models.
Display notation follows Relative Entropy Mixture.tex; CSV keys are internal.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import LogFormatterSciNotation, LogLocator, NullFormatter
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "report/data"
FIGURES = ROOT / "manuscript_theory/figs"
PAIRS = {
    "ising": (("I", "sigma"), ("sigma", "I"), ("I", "epsilon"),
              ("epsilon", "I"), ("sigma", "epsilon"), ("epsilon", "sigma")),
    "boson": (("00", "ss"), ("ss", "00"), ("00", "3s_s"),
              ("3s_s", "00"), ("ss", "3s_s"), ("3s_s", "ss")),
}
LABEL = {"I": "0", "sigma": r"\sigma", "epsilon": r"\varepsilon",
         "00": "00", "ss": "ss", "3s_s": "(3s,s)"}
STYLE = {"free": ("#0072B2", "-", "Free power"),
         "F1": ("#D55E00", "--", r"$M_1$"),
         "F2": ("#009E73", "-.", r"$M_2$"),
         "CFT": ("black", ":", "CFT")}
DATA_COLOR = "#666666"
MARKER = "o"
EXPECTED = {("small", "ising"): (0.02, 1255),
            ("small", "boson"): (0.1, 598),
            ("large", "ising"): (0.002, 260),
            ("large", "boson"): (0.005, 317)}


def read(name):
    with (DATA / name).open() as stream:
        return list(csv.DictReader(stream))


def predict(x, row):
    return sum(c * np.asarray(x) ** p for c, p in
               zip(json.loads(row["coefficients"]), json.loads(row["powers"])))


def key(row):
    return (int(row["N"]), int(row["ell"]))


def render(audit_path):
    plt.rcParams.update({"font.family": "serif", "font.size": 9,
        "mathtext.fontset": "cm", "axes.linewidth": 0.65,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True, "pdf.fonttype": 42,
        "savefig.bbox": "tight"})
    input_names = ["scaling_observations.csv", "scaling_fit_summary.csv",
                   "scaling_predictions.csv", "f2_summary_fits.csv"]
    hashes = {name: hashlib.sha256((DATA / name).read_bytes()).hexdigest()
              for name in input_names}
    observations = read(input_names[0])
    summaries = read(input_names[1])
    predictions = read(input_names[2])
    table_fits = {(r["region"], r["model"], r["a"], r["b"], r["fit"]): r
                  for r in read(input_names[3])}
    summary = {(r["region"], r["model"], r["a"], r["b"], r["fit"]): r
               for r in summaries}
    records = []
    for (region, model), (cutoff, count) in EXPECTED.items():
        groups = {}
        for a, b in PAIRS[model]:
            pts = sorted([r for r in observations if
                (r["region"], r["model"], r["a"], r["b"], r["fit_selected"])
                == (region, model, a, b, "1")], key=lambda r: float(r["x"]))
            assert len(pts) == count
            assert all(6 <= int(r["ell"]) <= 10 and 0 < float(r["x"]) <= cutoff
                       for r in pts)
            groups[a, b] = pts
            for fit in STYLE:
                row = summary[region, model, a, b, fit]
                table = table_fits[region, model, a, b, fit]
                assert int(row["m"]) == count and float(row["fit_max"]) == cutoff
                for field in ("powers", "coefficients"):
                    assert np.allclose(json.loads(row[field]), json.loads(table[field]),
                                       rtol=2e-13, atol=0), (model, region, fit, field)
        xmin = min(float(r["x"]) for pts in groups.values() for r in pts)
        for kind in ("value", "residual"):
            fig, axes = plt.subplots(3, 2, figsize=(7.1, 7.65), layout="constrained")
            for index, (a, b) in enumerate(PAIRS[model]):
                ax = axes.flat[index]
                pts = groups[a, b]
                x = np.array([float(r["x"]) for r in pts])
                q = np.array([float(r["value"]) for r in pts])
                xx = np.geomspace(xmin, cutoff, 600)
                if kind == "value":
                    ax.scatter(x, q, s=9, marker=MARKER, c=DATA_COLOR,
                               linewidths=0, alpha=0.7, zorder=2)
                    for fit, (color, linestyle, _) in STYLE.items():
                        yy = predict(xx, summary[region, model, a, b, fit])
                        assert np.isfinite(yy).all() and (yy > 0).all()
                        ax.plot(xx, yy, color=color, ls=linestyle, lw=1.25, zorder=3)
                    ax.set_yscale("log")
                else:
                    for fit, (color, _, _) in STYLE.items():
                        rr = sorted([r for r in predictions if
                            (r["region"], r["model"], r["a"], r["b"],
                             r["fit"], r["fit_selected"])
                            == (region, model, a, b, fit, "1")],
                            key=lambda r: float(r["x"]))
                        assert len(rr) == count and {key(r) for r in rr} == {key(r) for r in pts}
                        px = np.array([float(r["x"]) for r in rr])
                        pq = np.array([float(r["observed"]) for r in rr])
                        residual = np.array([float(r["relative_percent"]) for r in rr])
                        expected = 100 * (pq - predict(px, summary[region, model, a, b, fit])) / pq
                        assert np.allclose(residual, expected, rtol=1e-9, atol=2e-10)
                        ax.scatter(px, residual, s=8, marker=MARKER, c=color,
                                   linewidths=0, alpha=0.72, zorder=3)
                    ax.axhline(0, color=".6", lw=.6, zorder=1)
                    ax.set_yscale("linear")
                    ax.ticklabel_format(axis="y", style="plain", useOffset=False)
                    ax.margins(y=.12)
                ax.set_xscale("log")
                ax.set_xlim(xmin * .96, cutoff)
                ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5)))
                ax.xaxis.set_major_formatter(LogFormatterSciNotation(labelOnlyBase=False))
                ax.xaxis.set_minor_formatter(NullFormatter())
                xlabel = r"$x=l/R$" if region == "small" else r"$\epsilon=1-x=l/R$"
                observable = (r"$S(\rho_A^{\Psi}\Vert\bar{\rho}_A)$"
                              if region == "small" else
                              r"$\log 2-S(\rho_B^{\Psi}\Vert\bar{\rho}_B)$")
                ax.set_xlabel(xlabel)
                ax.set_title("(" + chr(97 + index) + ")  $("
                             + LABEL[a] + "," + LABEL[b] + ")$", loc="left")
                ax.set_ylabel(r"$R_i$ (%)" if kind == "residual" else observable)
                ax.tick_params(labelsize=8)
                for artist in ax.collections:
                    offsets = np.asarray(artist.get_offsets())
                    assert offsets.shape[0] == count and (offsets[:, 0] <= cutoff).all()
                    assert artist.get_paths().__len__() == 1
                assert ax.get_xlim()[1] == cutoff
                records.append(dict(region=region, model=model, kind=kind, pair=[a, b],
                    observations=count, x_min=float(x.min()), x_max=float(x.max()),
                    curve_max=float(xx[-1]), axis_max=float(ax.get_xlim()[1]),
                    marker=MARKER, length_color_encoding=False, data_color=DATA_COLOR,
                    residual_colors="model-specific; identical across all lengths",
                    yscale=ax.get_yscale(), excluded_points=0,
                    xlabel=xlabel, ylabel=ax.get_ylabel(),
                    state_pair=[LABEL[a], LABEL[b]],
                    legend_labels=[item[2] for item in STYLE.values()]))
            handles = []
            if kind == "value":
                handles.append(Line2D([], [], color=DATA_COLOR, marker=MARKER,
                                      ls="none", markersize=4, label="Numerics"))
            for fit, (color, linestyle, label) in STYLE.items():
                if fit == "CFT":
                    label = "CFT series" if region == "small" else "CFT leading"
                handles.append(Line2D([], [], color=color, ls=linestyle if kind == "value" else "none",
                    marker=None if kind == "value" else MARKER, markersize=4, label=label))
            fig.legend(handles=handles, loc="outside upper center", ncol=3 if kind == "value" else 4,
                       frameon=False, columnspacing=1.6)
            name = f"{model}_{region}_fit_window_{kind}.pdf"
            fig.savefig(FIGURES / name)
            plt.close(fig)
            print(f"Wrote {name}: {count} observations/pair, 0 < fraction <= {cutoff}", flush=True)
            audit_path.write_text(json.dumps({"input_sha256": hashes, "panels": records}, indent=2) + "\n")
    assert all(hashlib.sha256((DATA / name).read_bytes()).hexdigest() == sha
               for name, sha in hashes.items())
    assert len(records) == 48


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, required=True)
    render(parser.parse_args().audit)
