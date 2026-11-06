from pathlib import Path
import pandas as pd
from matplotlib.lines import Line2D
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Parameters
OUT_DIR = Path(__file__).resolve().parent / "results"
FIGSIZE = (3.8, 6.0)
FONTSIZE = 11
LABEL_SIZE = 9
COLORS = {'neither': '#C5C9CE', 'DAS only': '#C75559', 'DEGs only': '#4379A6', 'both': '#9366AC'}
plt.rcParams.update({'font.family': 'Arial', 'font.size': FONTSIZE, 'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none', 'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic', 'axes.linewidth': 0.8, 'xtick.labelsize': 10, 'ytick.labelsize': 10})
OUT_DIR.mkdir(exist_ok=True)

# Load data
source = pd.read_csv(OUT_DIR / 'four_class_points.tsv.gz', sep='\t', low_memory=False)
source = source[source.plottable].copy()
# 四分类来自先完成的统计表
labels = pd.DataFrame([
    ('AD_vs_Young', 'Tmem214', 'TH', 2.32, 1.0),
    ('AD_vs_Young', 'Dnm1l', 'CTX_2', 2.32, 0.77),
    ('AD_vs_Young', 'Clta', 'FT', -2.32, 0.6),
    ('AD_vs_Young', 'Kbtbd3', 'HIP', 2.32, 0.48),
    ('AD_vs_Young', 'Mbp', 'CA1', -2.32, 0.23),
    ('AD_vs_Young', 'Cd74', 'CTX_3', 2.32, 0.26),
    ('Aged_vs_Young', 'Atp5j2', 'CA3', 2.32, 1.0),
    ('Aged_vs_Young', 'Ncam1', 'FT', -2.32, 0.85),
    ('Aged_vs_Young', 'Mag', 'CTX_2', 2.32, 0.62),
    ('Aged_vs_Young', 'Rps24', 'FT', -2.32, 0.43),
    ('Aged_vs_Young', 'Zcrb1', 'FT', 2.32, 0.39),
    ('Aged_vs_Young', 'Mbp', 'CA1', -2.32, 0.22),
    ('Aged_vs_Young', 'C4b', 'CTX_3', 2.32, 0.2),
], columns=["comparison", "gene_name", "region", "label_x", "label_y"])
labels = labels.merge(source, on=['comparison', 'gene_name', 'region'], how='left', validate='one_to_one')
assert labels.plot_y.notna().all()
assert source.n_testable.ge(1).all()
assert source.loc[source.panel_class.isin(['DEGs only', 'both']), 'padj'].lt(0.05).all()
summaries = []

# Plot
for layer, prefix in [('AS', 'AS_gene_region_DEG_DAS_four_classes')]:
    name = 'AS'
    ylabel = 'Max |ΔPSI|'
    data = source.copy()
    data.to_csv(OUT_DIR / f'{prefix}_source.tsv.gz', sep='\t', index=False)
    lab_all = labels.copy()
    lab_all.to_csv(OUT_DIR / f'{prefix}_labels.tsv', sep='\t', index=False)
    fig, axes = plt.subplots(2, 1, figsize=FIGSIZE)
    fig.subplots_adjust(left=0.19, right=0.975, bottom=0.24, top=0.935, hspace=0.43)
    for row, comparison in enumerate(['AD_vs_Young', 'Aged_vs_Young']):
        ax = axes[row]
        df = data[data.comparison.eq(comparison)]
        for cls in COLORS:
            sub = df[df.panel_class.eq(cls)]
            ax.scatter(sub.log2FC, sub.plot_y, s=5 if cls == 'neither' else 7, color=COLORS[cls], alpha=0.45 if cls == 'neither' else 0.55, linewidths=0, edgecolors='none', rasterized=True, zorder=2 if cls == 'neither' else 3)
            summaries.append({'comparison': comparison, 'layer': layer, 'category': cls, 'n_region_gene': len(sub), 'n_in_focus': int(sub.in_focus.sum())})
        ax.axvline(-0.5, color='black', linestyle='--', linewidth=0.8, zorder=1)
        ax.axvline(0.5, color='black', linestyle='--', linewidth=0.8, zorder=1)
        ax.axhline(0.05, color='black', linestyle='--', linewidth=0.8, zorder=1)
        ax.set(xlim=(-2.5, 2.5), ylim=(-0.02, 1.09), xlabel='Gene log2FC', ylabel=ylabel)
        ax.set_xticks([-2, -1, 0, 1, 2])
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
        ax.set_title(('AD vs Young' if row == 0 else 'Aged vs Young') + ' · ' + name, fontsize=11, pad=10)
        ax.text(-0.22, 1.055, 'ab'[row], transform=ax.transAxes, fontsize=12, weight='bold')
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(length=3, pad=3)
        ax.grid(False)
        lab = lab_all[lab_all.comparison.eq(comparison)]
        for _, r in lab.iterrows():
            ax.annotate(r.gene_name + ' | ' + r.region, xy=(r.log2FC, r.plot_y), xytext=(-2.32 if r.label_x < 0 else 2.32, r.label_y), ha='left' if r.label_x < 0 else 'right', va='center', fontsize=LABEL_SIZE, color='#252525', arrowprops={'arrowstyle': '-', 'color': '#686868', 'lw': 0.55, 'shrinkA': 2, 'shrinkB': 3}, zorder=6)
    legend_labels = ['DAS only', 'DEGs only', 'both', 'neither']
    handles = [Line2D([], [], marker='o', linestyle='none', color=COLORS[k], markeredgewidth=0, markersize=5, label=t) for k, t in zip(legend_labels, legend_labels)]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.56, 0.025), frameon=False, fontsize=9, handletextpad=0.3, labelspacing=0.5)
    fig.savefig(OUT_DIR / f'{prefix}.pdf', dpi=600)
    fig.savefig(OUT_DIR / f'{prefix}.svg', dpi=600)
    fig.savefig(OUT_DIR / f'{prefix}.png', dpi=300)
    # Vector output
    for ax in axes:
        for collection in ax.collections:
            collection.set_rasterized(False)
    fig.savefig(OUT_DIR / f'{prefix}_vector.pdf')
    # RGB散点层：去除透明蒙版，保留矢量文字、坐标轴和引线
    import numpy as np
    for row, ax in enumerate(axes):
        box = ax.get_position()
        temp, ta = plt.subplots(figsize=(FIGSIZE[0]*box.width, FIGSIZE[1]*box.height), dpi=600)
        temp.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ta.set(xlim=ax.get_xlim(), ylim=ax.get_ylim())
        ta.set_axis_off()
        df = data[data.comparison.eq(['AD_vs_Young', 'Aged_vs_Young'][row])]
        for cls in COLORS:
            sub = df[df.panel_class.eq(cls)]
            ta.scatter(sub.log2FC, sub.plot_y, s=5 if cls == 'neither' else 7, color=COLORS[cls], alpha=0.45 if cls == 'neither' else 0.55, linewidths=0, edgecolors='none')
        temp.canvas.draw()
        rgb = np.asarray(temp.canvas.buffer_rgba())[:, :, :3].copy()
        plt.close(temp)
        for collection in list(ax.collections):
            collection.remove()
        for patch in list(ax.patches):
            patch.remove()
        ax.imshow(rgb, extent=(*ax.get_xlim(), *ax.get_ylim()), origin='upper', aspect='auto', interpolation='none', zorder=0)
    fig.savefig(OUT_DIR / f'{prefix}_AI_compatible.pdf', dpi=600)
    fig.savefig(OUT_DIR / f'{prefix}_AI_compatible.png', dpi=300)
    fig.savefig(OUT_DIR / f'{prefix}_AI_compatible.svg', dpi=600)
    plt.close(fig)

summary = pd.DataFrame(summaries)
summary.to_csv(OUT_DIR / 'category_counts.tsv', sep='\t', index=False)
print(summary.to_string(index=False))
expected = pd.read_csv(OUT_DIR / 'four_class_counts.tsv', sep='\t')
expected = expected[expected.region.eq('ALL')]
check = summary.merge(expected, on=['comparison', 'category'], suffixes=('_plot', '_stats'), validate='one_to_one')
assert len(check) == 8
assert (check.n_region_gene_plot == check.n_region_gene_stats).all()
assert (check.n_in_focus_plot == check.n_in_focus_stats).all()
print('Verified: all four plot categories match statistics, including visible points.')
