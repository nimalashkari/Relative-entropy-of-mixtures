"""Render all equal-mixture panels from complete step-4 observations and saved fits."""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/entropy-step4-matplotlib')
os.environ.setdefault('XDG_CACHE_HOME','/tmp/entropy-step4-cache')
import csv
import hashlib
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import LogFormatterSciNotation, LogLocator, NullFormatter
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
MANUSCRIPT = ROOT/'manuscript_theory'
DATA = ROOT/'cluster_uniform_step4/results/full_observations.csv'
FITS = MANUSCRIPT/'numerical/step4/fits.csv'
OUT = MANUSCRIPT/'figs'
AUDIT = MANUSCRIPT/'numerical/step4/figure_audit.json'
EXPECTED = {('small','ising'):(.02,9744,1974),('small','boson'):(.02,19984,4022),
            ('large','ising'):(.005,8245,1749),('large','boson'):(.005,18485,3797)}
LABEL = {'I':'0','sigma':r'\sigma','epsilon':r'\varepsilon','00':'00','ss':'ss','3s_s':'(3s,s)'}
STYLE = {'free':('#0072B2','-','Free power'), 'M1':('#D55E00','--',r'$M_1$'),
         'M2':('#009E73','-.',r'$M_2$'), 'CFT':('black',':','CFT (known terms)')}


def predict(x,r,kind):
    if kind == 'free':
        return float(r['C_free'])*x**float(r['r_free'])
    c0 = float(r['c0_cft'] if kind == 'CFT' else r['c0_m1'] if kind == 'M1' else r['c0_m2'])
    value = c0*x**float(r['r0'])
    if kind == 'M2' or (kind == 'CFT' and r['c1_cft']):
        value += float(r['c1_cft'] if kind == 'CFT' else r['c1_m2'])*x**float(r['r1'])
    return value


def main():
    digest = hashlib.sha256(DATA.read_bytes()).hexdigest()
    assert digest == '32c3a3d54eee116a7d5e17d48d6074fa92c41d6bcf7e59dbb6eb889412ead5e9'
    with FITS.open() as f:
        fits = list(csv.DictReader(f))
    grouped = {(r['region'],r['model'],r['a'],r['b']):[] for r in fits}
    with DATA.open() as f:
        for r in csv.DictReader(f):
            key = tuple(r[k] for k in ('region','model','a','b'))
            if float(r['x']) <= EXPECTED[key[:2]][0]:
                grouped[key].append(r)
    plt.rcParams.update({'font.family':'serif','font.size':9,'mathtext.fontset':'cm',
                         'axes.linewidth':.65,'xtick.direction':'in','ytick.direction':'in',
                         'xtick.top':True,'ytick.right':True,'pdf.fonttype':42,'savefig.bbox':'tight'})
    records = []
    for (region,model),(cap,count,sizes) in EXPECTED.items():
        local = [r for r in fits if (r['region'],r['model']) == (region,model)]
        xmin = min(float(p['x']) for r in local for p in grouped[region,model,r['a'],r['b']])
        for kind in ('value','residual'):
            fig,axes = plt.subplots(3,2,figsize=(7.1,7.65),layout='constrained')
            for index,r in enumerate(local):
                ax = axes.flat[index]
                points = sorted(grouped[region,model,r['a'],r['b']],key=lambda p:(float(p['x']),int(p['N']),int(p['ell'])))
                assert len(points) == int(r['n']) == count
                assert len({p['N'] for p in points}) == int(r['sizes']) == sizes
                assert all(int(p['N'])%4==0 and 6<=int(p['ell'])<=10 for p in points)
                x = np.array([float(p['x']) for p in points]); y = np.array([float(p['value']) for p in points])
                xx = np.geomspace(xmin,cap,600)
                residual_ranges = {}
                if kind == 'value':
                    ax.scatter(x,y,s=2,marker='o',c='#666666',linewidths=0,alpha=.35,rasterized=True,zorder=2)
                for fit,(color,style,_) in STYLE.items():
                    yy = predict(xx,r,fit)
                    assert np.all(np.isfinite(yy)) and np.all(yy>0)
                    residual = 100*(y-predict(x,r,fit))/y
                    assert np.isfinite(residual).all()
                    residual_ranges[fit] = [float(residual.min()),float(residual.max())]
                    if kind == 'value':
                        ax.plot(xx,yy,color=color,ls=style,lw=1.25,zorder=3)
                    else:
                        ax.scatter(x,residual,s=2,marker='o',c=color,linewidths=0,alpha=.4,rasterized=True,zorder=3)
                if kind == 'value':
                    ax.set_yscale('log')
                else:
                    ax.axhline(0,color='.6',lw=.6,zorder=1)
                    ax.ticklabel_format(axis='y',style='plain',useOffset=False)
                    ax.margins(y=.12)
                ax.set_xscale('log'); ax.set_xlim(xmin*.96,cap)
                ax.xaxis.set_major_locator(LogLocator(base=10,subs=(1,2,5)))
                ax.xaxis.set_major_formatter(LogFormatterSciNotation(labelOnlyBase=False))
                ax.xaxis.set_minor_formatter(NullFormatter())
                ax.set_xlabel(r'$x=l/R$' if region=='small' else r'$\epsilon=1-x=l/R$')
                ax.set_ylabel(r'$R_i$ (%)' if kind=='residual' else r'$S(\rho_A^\Psi\Vert\bar\rho_A)$' if region=='small' else r'$\log 2-S(\rho_B^\Psi\Vert\bar\rho_B)$')
                ax.set_title('('+chr(97+index)+') $('+LABEL[r['a']]+','+LABEL[r['b']]+')$',loc='left')
                ax.tick_params(labelsize=8)
                assert all(len(artist.get_offsets())==count for artist in ax.collections)
                records.append(dict(region=region,model=model,a=r['a'],b=r['b'],kind=kind,
                                    observations=count,sizes=sizes,minimum_fraction=float(x.min()),maximum_fraction=float(x.max()),
                                    upper_cutoff=cap,lengths=sorted({int(p['ell']) for p in points}),
                                    cft_has_next_term=bool(r['c1_cft']),residual_range_percent=residual_ranges,
                                    plot_data_sha256=hashlib.sha256(np.column_stack([x,y]).tobytes()).hexdigest()))
            handles = []
            if kind == 'value':
                handles.append(Line2D([],[],color='#666666',marker='o',ls='none',markersize=4,label='Numerics'))
            handles += [Line2D([],[],color=color,ls=style if kind=='value' else 'none',marker=None if kind=='value' else 'o',markersize=4,label=name)
                        for color,style,name in STYLE.values()]
            fig.legend(handles=handles,loc='outside upper center',ncol=3 if kind=='value' else 4,frameon=False,columnspacing=1.4)
            name = f'{model}_{region}_step4_{kind}.pdf'
            fig.savefig(OUT/name,dpi=350)
            plt.close(fig)
            AUDIT.write_text(json.dumps(dict(source_sha256=digest,fits_sha256=hashlib.sha256(FITS.read_bytes()).hexdigest(),panels=records),indent=2)+'\n')
            print(f'{name}: {count} observations per direction; all points rendered',flush=True)
    assert len(records) == 48


if __name__ == '__main__':
    main()
