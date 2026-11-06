from pathlib import Path
import warnings
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
import scanpy as sc

# Parameters
ROOT = Path(__file__).resolve().parents[1]
SCORE_DIR = ROOT / 'reference_results'  # Change to outputs/03_scores after recalculation.
OUT = ROOT / 'outputs/figures'
SELECTED = ['G21-4', 'G21-8', 'G21-5']
SAMPLES = ['G21-4', 'G21-1', 'G21-2', 'G21-7', 'G21-8', 'G21-5']
GROUP = {'G21-4': 'Young', 'G21-7': 'Young', 'G21-1': 'Aged', 'G21-8': 'Aged', 'G21-2': 'AD', 'G21-5': 'AD'}
REGIONS = ['CTX_1', 'CTX_2', 'CTX_3', 'HIP', 'CA1', 'CA3', 'DG', 'TH', 'HY', 'STR', 'FT', 'VS', 'EPI']
LIMITS = {'Aging': (-0.12, 0.22), 'AD': (-0.13, 0.33)}
COLORS = {'Aging': ['#527D9A', '#F7F5F0', '#B74626'], 'AD': ['#478B83', '#F7F5F0', '#74418B']}
SPOT_SIZE = 25
DPI = 300
plt.rcParams.update({'font.family': 'Arial', 'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none', 'font.size': 8, 'axes.titlesize': 9, 'axes.labelsize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7, 'axes.grid': False})


def spatial(ax, frame, column, norm, title, cmap):
    a = ad.AnnData(X=sparse.csr_matrix((len(frame), 1)), obs=frame.set_index('spot_id').copy())
    a.obsm['spatial'] = frame[['x', 'y']].to_numpy()
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='Use `squidpy.pl.spatial_scatter` instead.', category=FutureWarning)
        sc.pl.spatial(a, color=column, ax=ax, spot_size=SPOT_SIZE, img_key=None, norm=norm, cmap=cmap, colorbar_loc=None, frameon=False, show=False, title=title, na_color='#BDBDBD', na_in_legend=False)
    ax.set_xlabel(''); ax.set_ylabel(''); ax.grid(False)


def save(fig, name):
    fig.savefig(OUT / f'{name}.pdf', bbox_inches='tight')
    fig.savefig(OUT / f'{name}.svg', bbox_inches='tight')
    fig.savefig(OUT / f'{name}.png', dpi=DPI, bbox_inches='tight')
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    print('Plot score source:', SCORE_DIR)
    spots = pd.read_csv(SCORE_DIR / 'spot_scores.tsv.gz', sep='\t')
    regional = pd.read_csv(SCORE_DIR / 'sample_region_scores.csv')
    regional = regional.loc[regional.region.isin(REGIONS)].copy()
    d = spots[['sample', 'spot_id', 'Regions', 'x', 'y']].copy()
    for score in ['Aging', 'AD']:
        index = regional.loc[regional.score.eq(score)].set_index(['sample', 'region']).adjusted_mean
        d[score + '_region'] = pd.MultiIndex.from_frame(d[['sample', 'Regions']]).map(index).to_numpy()
    regional.to_csv(OUT / 'region_means_source.csv', index=False)
    d.to_csv(OUT / 'regional_plot_values.tsv.gz', sep='\t', index=False)
    pd.DataFrame([{'score': s, 'vmin': lo, 'vmax': hi} for s, (lo, hi) in LIMITS.items()]).to_csv(OUT / 'color_limits.csv', index=False)
    cmaps = {s: LinearSegmentedColormap.from_list(s, c) for s, c in COLORS.items()}
    fig = plt.figure(figsize=(7.2047244, 5.2))
    gs = fig.add_gridspec(2, 4, width_ratios=[1, 1, 1, .055], left=.03, right=.91, bottom=.085, top=.94, wspace=.10, hspace=.15)
    for row, score in enumerate(['Aging', 'AD']):
        norm = Normalize(*LIMITS[score])
        for j, sample in enumerate(SELECTED):
            ax = fig.add_subplot(gs[row, j])
            spatial(ax, d.loc[d['sample'].eq(sample)], score + '_region', norm, f'{GROUP[sample]} · {sample}' if row == 0 else '', cmaps[score])
        cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmaps[score]), cax=fig.add_subplot(gs[row, 3]), extend='neither')
        cb.set_label(score + ' score − Young baseline\nRegion mean')
    fig.text(.03, .038, 'Region means of measured spots; within-region variation is not shown.', fontsize=6.5)
    fig.text(.03, .018, 'Grey: unavailable baseline or region summary.', fontsize=6.5)
    save(fig, 'Fig2_region_means')
    # Retain all biological replicates in supplementary spatial maps.
    for score in ['Aging', 'AD']:
        fig, axes = plt.subplots(2, 3, figsize=(7.2047244, 5.2))
        fig.subplots_adjust(left=.02, right=.875, bottom=.07, top=.92, wspace=.04, hspace=.13)
        norm = Normalize(*LIMITS[score])
        for ax, sample in zip(axes.flat, SAMPLES):
            spatial(ax, d.loc[d['sample'].eq(sample)], score + '_region', norm, f'{GROUP[sample]} · {sample}', cmaps[score])
        cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmaps[score]), cax=fig.add_axes([.90, .25, .018, .43]), extend='neither')
        cb.set_label(score + ' score − Young baseline\nRegion mean')
        fig.text(.02, .025, 'Region means of measured spots.', fontsize=6.5)
        save(fig, f'FigS_all_six_{score}_region_means')
    print(regional[['score', 'sample', 'region', 'n_spots', 'adjusted_mean']].to_string(index=False))
    print('Measured spots:', len(d))


if __name__ == '__main__':
    main()
