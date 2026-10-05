"""Manuscript figures from frozen M2 coefficients; no numerical refitting."""
from pathlib import Path
import csv
import json
import hashlib
import argparse
import numpy as np
from scipy.special import loggamma, digamma, polygamma
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

MANUSCRIPT = Path(__file__).resolve().parents[2]
DATA = MANUSCRIPT / 'numerical/y_step4'
FIG = MANUSCRIPT / 'figs'
MARKER_ALPHA = .45
RESIDUAL_ALPHA = .55
LABEL = {'I': '0', 'sigma': r'\sigma', 'epsilon': r'\varepsilon',
         '00': '00', 'ss': 'ss', '3s_s': '(3s,s)'}
PAIRS = {'boson': [('00', 'ss'), ('00', '3s_s'), ('ss', '3s_s')],
         'ising': [('I', 'sigma'), ('I', 'epsilon'), ('sigma', 'epsilon')]}

def mu(delta, y):
    z = delta + 1 + 1j * np.log((1-y)/y) / (2*np.pi)
    return np.sqrt((1-y)/y) * np.exp(2*loggamma(z).real-loggamma(2*delta+2)) * (.5+digamma(z).imag/np.pi)

def read(name):
    return list(csv.DictReader((DATA/name).open()))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--preview-dir',type=Path,default=DATA/'preview')
    args=parser.parse_args();args.preview_dir.mkdir(parents=True,exist_ok=True)
    rows = read('amplitude_comparison.csv')
    fits = read('fits.csv')
    plt.rcParams.update({'font.size': 11, 'axes.titlesize': 12,
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.grid': True, 'grid.alpha': .17, 'pdf.fonttype': 42,
        'savefig.facecolor': 'white'})
    audit = {'data_sha256': hashlib.sha256((DATA/'amplitude_comparison.csv').read_bytes()).hexdigest(), 'groups': []}
    for model in ('boson', 'ising'):
        for region in ('small', 'large'):
            fig, axes = plt.subplots(2, 3, figsize=(8.8, 4.5), layout='constrained',
                                     gridspec_kw={'height_ratios': [2, 1]})
            yy = np.linspace(.4, .6, 401)
            for j, (a, b) in enumerate(PAIRS[model]):
                ax, rx = axes[:, j]
                center = next(r for r in fits if (r['model'], r['region'], r['a'], r['b'], float(r['y']), r['fit']) == (model, region, a, b, .5, 'M2'))
                half = float(center['c0_cft']); delta = float(center['r0'])/2
                cft = half * (4*(1-yy)**2 if region == 'small' else mu(delta, yy)/mu(delta, .5))
                ax.plot(yy, cft, color='#252b32', lw=1.6)
                if region == 'large':
                    aa = 2 + 4*polygamma(1, delta+1)/np.pi**2
                    ax.plot(yy, half*(1-aa*(yy-.5)+aa*(yy-.5)**2), color='#bd8a26', ls='--', lw=1.5)
                for reverse, (u, v) in enumerate(((a, b), (b, a))):
                    data = sorted((r for r in rows if (r['model'], r['region'], r['a'], r['b']) == (model, region, u, v)), key=lambda r: float(r['y']))
                    y = np.array([float(r['y']) for r in data])
                    value = np.array([float(r['c0_fit']) for r in data])
                    theory = np.array([float(r['c0_theory']) for r in data])
                    error = 100*(value/theory-1)
                    assert len(data) == 40 and np.allclose(y, np.linspace(.4, .6, 40), rtol=0, atol=2e-16)
                    predicted = half*(4*(1-y)**2 if region == 'small' else mu(delta, y)/mu(delta, .5))
                    assert np.allclose(theory, predicted, rtol=2e-14)
                    for r, e in zip(data, error):
                        fit = next(f for f in fits if (f['model'], f['region'], f['a'], f['b'], f['y'], f['fit']) == (model, region, u, v, r['y'], 'M2'))
                        assert float(fit['c0']) == float(r['c0_fit'])
                        assert abs(e-float(fit['c0_error_pct'])) < 1e-10
                    color, marker = ('#2369a8', 'o') if not reverse else ('#c05a24', 's')
                    ax.scatter(y, value, s=13, marker=marker, facecolors=color if not reverse else 'none', edgecolors=color, lw=.8, alpha=MARKER_ALPHA, zorder=4)
                    rx.plot(y, error, color=color, marker=marker, ms=2.4, lw=.8, mfc=color if not reverse else 'none', alpha=RESIDUAL_ALPHA)
                    audit['groups'].append(dict(model=model, region=region, a=u, b=v, points=len(data), max_abs_eta0=float(max(abs(error)))))
                ax.set_title('$'+LABEL[a]+r'\ \leftrightarrow\ '+LABEL[b]+'$')
                ax.set_ylabel(r'$\bar c_0(y),\ c_0^{\rm CFT}(y)$' if j == 0 else '')
                rx.set_ylabel(r'$\eta_0\ [\%]$' if j == 0 else '')
                rx.set_xlabel('$y$'); rx.axhline(0, color='black', lw=.65)
                for panel in (ax, rx):
                    panel.set_xlim(.395, .605)
                    panel.set_xticks([.4, .45, .5, .55, .6])
                    panel.tick_params(labelsize=9)
                    panel.ticklabel_format(axis='y', style='plain', useOffset=False)
                ax.tick_params(labelbottom=False)
            handles = [Line2D([], [], color='#252b32', lw=1.6, label=r'$c_0^{\rm CFT}(y)$'),
                       Line2D([], [], color='#2369a8', marker='o', ls='', ms=4, alpha=MARKER_ALPHA, label=r'$\bar c_0(y)$: left to right'),
                       Line2D([], [], color='#c05a24', marker='s', mfc='none', ls='', ms=4, alpha=MARKER_ALPHA, label=r'$\bar c_0(y)$: right to left')]
            if region == 'large':
                handles.append(Line2D([], [], color='#bd8a26', ls='--', lw=1.5, label='Quadratic expansion'))
            fig.legend(handles=handles, loc='outside lower center', ncol=2, fontsize=10, frameon=False)
            name = f'{model}_{region}_general_y_step4_coefficients'
            fig.savefig(FIG/(name+'.pdf'), bbox_inches='tight')
            fig.savefig(args.preview_dir/(name+'.png'), dpi=180, bbox_inches='tight')
            plt.close(fig)
    assert sum(g['points'] for g in audit['groups']) == 960
    (DATA/'figure_verification.json').write_text(json.dumps(audit, indent=2)+'\n')
    print('Verified and plotted all 960 coefficients and eta0 values in four figures.')

if __name__ == '__main__':
    main()
