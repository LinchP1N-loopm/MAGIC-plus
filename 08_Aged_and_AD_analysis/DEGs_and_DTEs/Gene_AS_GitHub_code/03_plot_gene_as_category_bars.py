from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Parameters
OUT = Path(__file__).resolve().parent / "results"
ORDER = ['DAS only', 'both', 'DEGs only', 'neither']
COLORS = ['#C75559', '#9366AC', '#4379A6', '#C5C9CE']
FIGSIZE = (4.2, 3.6)
FONTSIZE = 11
plt.rcParams.update({'font.family': 'Arial', 'font.size': FONTSIZE, 'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none', 'axes.linewidth': 1.1})

# Load data
counts = pd.read_csv(OUT / 'four_class_counts.tsv', sep='\t')
data = counts[counts.region.eq('ALL')].copy()
data['percentage'] = data.n_region_gene / data.groupby('comparison').n_region_gene.transform('sum') * 100
data.to_csv(OUT / 'category_bars_source.tsv', sep='\t', index=False)

# Plot
combined, axes = plt.subplots(1, 2, figsize=(8.4, 3.6))
for i, (comparison, title) in enumerate([('AD_vs_Young', 'AD vs Young'), ('Aged_vs_Young', 'Aged vs Young')]):
    d = data[data.comparison.eq(comparison)].set_index('category').loc[ORDER]
    fig, ax = plt.subplots(figsize=FIGSIZE)
    for target in [ax, axes[i]]:
        target.bar(range(4), d.n_region_gene, width=0.72, color=COLORS, edgecolor='#444444', linewidth=0.7)
        for x, (_, row) in enumerate(d.iterrows()):
            target.text(x, row.n_region_gene + 450, f'{int(row.n_region_gene):,}\n({row.percentage:.1f}%)', ha='center', va='bottom', fontsize=10)
        target.set_xticks(range(4), ORDER)
        target.set(xlim=(-0.55, 3.55), ylim=(0, 35000), ylabel='Number of\n(gene, brain region) pairs', xlabel='Test', title=title)
        target.set_yticks([0, 5000, 10000, 15000, 20000, 25000, 30000])
        target.tick_params(width=1, length=3)
        target.grid(False)
    fig.tight_layout(pad=0.7)
    fig.savefig(OUT / f'category_bars_{comparison}.pdf')
    fig.savefig(OUT / f'category_bars_{comparison}.svg')
    fig.savefig(OUT / f'category_bars_{comparison}.png', dpi=600)
    plt.close(fig)
combined.tight_layout(pad=0.7, w_pad=1.7)
combined.savefig(OUT / 'category_bars_combined.pdf')
combined.savefig(OUT / 'category_bars_combined.svg')
combined.savefig(OUT / 'category_bars_combined.png', dpi=600)
plt.close(combined)
assert data.groupby('comparison').percentage.sum().sub(100).abs().lt(1e-8).all()
print(data[['comparison', 'category', 'n_region_gene', 'percentage']].to_string(index=False))
