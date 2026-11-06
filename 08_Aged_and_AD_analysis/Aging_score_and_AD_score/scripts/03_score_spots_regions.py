from pathlib import Path
import json
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr

# Parameters
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data/h5ad'
SIGNATURES = ROOT / 'data/fixed_signatures'  # Use outputs/02_signatures to score rebuilt signatures.
OUT = ROOT / 'outputs/03_scores'
GROUP_MAP = {'G21-4': 'Young', 'G21-7': 'Young', 'G21-1': 'Aged', 'G21-8': 'Aged', 'G21-2': 'AD', 'G21-5': 'AD'}
REGIONS = ['CTX_1', 'CTX_2', 'CTX_3', 'HIP', 'CA1', 'CA3', 'DG', 'TH', 'HY', 'STR', 'FT', 'VS']

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    var = pd.read_csv(SIGNATURES / 'common_genes.csv').set_index('gene_id')
    signatures = {name: pd.read_csv(SIGNATURES / f'{name}_signature.csv') for name in ['Aging', 'AD']}
    controls = pd.read_csv(SIGNATURES / 'module_controls.csv')
    genes = var.index.to_numpy()
    modules, status = {}, []
    for name, sig in signatures.items():
        sizes = sig.direction.value_counts()
        status.append({'score': name, 'available': True, 'Up': int(sizes['Up']), 'Down': int(sizes['Down'])})
        for direction in ['Up', 'Down']:
            target = sig.loc[sig.direction.eq(direction), 'gene_id']
            background = controls.loc[controls.module.eq(f'{name}_{direction}'), 'gene_id']
            idx, ctrl = var.index.get_indexer(target), var.index.get_indexer(background)
            assert len(idx) and len(ctrl) and (idx >= 0).all() and (ctrl >= 0).all()
            modules[f'{name}_{direction}'] = (idx, ctrl)
    frames = []
    for sample, condition in GROUP_MAP.items():
        a = ad.read_h5ad(DATA / f'{sample}_regions.h5ad', backed='r')
        o = a.obs[['Regions']].copy()
        o['sample'], o['condition'], o['spot_id'] = sample, condition, a.obs_names.astype(str)
        o['x'], o['y'] = a.obsm['spatial'][:, 0], a.obsm['spatial'][:, 1]
        counts = sparse.csr_matrix(a.layers['counts'])
        o['total_counts'] = np.asarray(counts.sum(axis=1)).ravel()
        o['n_genes'] = np.asarray((counts > 0).sum(axis=1)).ravel()
        frames.append(o.reset_index(drop=True))
        a.file.close()
    spots = pd.concat(frames, ignore_index=True)
    detection, contribs = [], []
    for sample in GROUP_MAP:
        a = ad.read_h5ad(DATA / f'{sample}_regions.h5ad', backed='r')
        feature_idx = a.var_names.get_indexer(genes)
        assert (feature_idx >= 0).all(), 'Missing scoring feature; use the original feature space.'
        x = sparse.csr_matrix(a.X[:, feature_idx])
        mask = spots['sample'].eq(sample)
        assert spots.loc[mask, 'spot_id'].tolist() == a.obs_names.tolist()
        for key, (idx, ctrl) in modules.items():
            target = np.asarray(x[:, idx].mean(axis=1)).ravel()
            bg = np.asarray(x[:, ctrl].mean(axis=1)).ravel()
            spots.loc[mask, key] = target - bg
            spots.loc[mask, f'{key}_detected'] = np.asarray((x[:, idx] > 0).sum(axis=1)).ravel()
            detection.append({'sample': sample, 'module': key, 'n_genes': len(idx), 'mean_detected': float(np.asarray((x[:, idx] > 0).sum(axis=1)).mean()), 'fraction_zero_detected': float(np.mean(np.asarray((x[:, idx] > 0).sum(axis=1)).ravel() == 0))})
            contribs.extend({'sample': sample, 'condition': GROUP_MAP[sample], 'module': key, 'gene_id': genes[i], 'symbol': var.iloc[i].Symbol, 'mean_X': float(v)} for i, v in zip(idx, np.asarray(x[:, idx].mean(axis=0)).ravel()))
        a.file.close()
        del a, x
        print('Scored', sample, flush=True)
    pd.DataFrame(detection).to_csv(OUT / 'module_detection.csv', index=False)
    pd.DataFrame(contribs).to_csv(OUT / 'signature_gene_sample_expression.csv', index=False)
    baseline_rows, effects, correlations = [], [], []
    available_scores = [s['score'] for s in status if s['available']]
    for name in available_scores:
        col = f'{name}Score'
        spots[col] = spots[f'{name}_Up'] - spots[f'{name}_Down']
        # 指定12脑区都有两只Young；额外区域仅在两只Young均有时校正
        young = spots.loc[spots.condition.eq('Young')].groupby(['sample', 'Regions'])[col].median().reset_index()
        b = young.groupby('Regions')[col].agg(['mean', 'count'])
        b.loc[b['count'] != 2, 'mean'] = np.nan
        spots[f'{col}_adjusted'] = spots[col] - spots.Regions.map(b['mean'])
        baseline_rows.extend({'score': name, 'region': r, 'Young_baseline': row['mean'], 'n_Young': int(row['count'])} for r, row in b.iterrows())
        for sample, d in spots.groupby('sample'):
            for cov in ['total_counts', 'n_genes']:
                correlations.append({'score': name, 'sample': sample, 'covariate': cov, 'spearman_r': spearmanr(d[col], d[cov]).statistic})
    pd.DataFrame(baseline_rows).to_csv(OUT / 'Young_region_baselines.csv', index=False)
    pd.DataFrame(correlations).to_csv(OUT / 'score_QC_correlations.csv', index=False)
    spots.to_csv(OUT / 'spot_scores.tsv.gz', sep='\t', index=False)
    summaries = []
    for name in available_scores:
        col = f'{name}Score'
        for (sample, condition, region), d in spots.groupby(['sample', 'condition', 'Regions']):
            summaries.append({'score': name, 'sample': sample, 'condition': condition, 'region': region, 'n_spots': len(d), 'mean': d[col].mean(), 'median': d[col].median(), 'adjusted_mean': d[f'{col}_adjusted'].mean(), 'adjusted_median': d[f'{col}_adjusted'].median()})
    summary = pd.DataFrame(summaries)
    summary.to_csv(OUT / 'sample_region_scores.csv', index=False)
    whole = []
    for name in available_scores:
        for (sample, condition), d in spots.groupby(['sample', 'condition']):
            whole.append({'score': name, 'sample': sample, 'condition': condition, 'n_spots': len(d), 'mean': d[f'{name}Score'].mean(), 'median': d[f'{name}Score'].median(), 'region_equal_mean': summary.loc[summary.score.eq(name) & summary['sample'].eq(sample) & summary.region.isin(REGIONS), 'mean'].mean()})
        for region in REGIONS:
            r = summary.loc[summary.score.eq(name) & summary.region.eq(region)]
            y = r.loc[r.condition.eq('Young'), 'mean'].to_numpy()
            for condition in ['Aged', 'AD']:
                t = r.loc[r.condition.eq(condition), 'mean'].to_numpy()
                effects.append({'score': name, 'region': region, 'comparison': f'{condition}_vs_Young', 'difference_of_mouse_means': float(t.mean() - y.mean()), 'n_positive_pairwise': int((t[:, None] > y).sum()), 'n_pairwise': 4})
    pd.DataFrame(whole).to_csv(OUT / 'sample_scores.csv', index=False)
    pd.DataFrame(effects).to_csv(OUT / 'region_effects.csv', index=False)
    (OUT / 'parameters.json').write_text(json.dumps({'signatures': str(SIGNATURES.relative_to(ROOT)), 'input': 'existing normalized log1p X; no renormalization', 'controls': 'loaded fixed module_controls.csv', 'baseline': 'mean of two Young sample-region medians'}, indent=2))
    print(pd.DataFrame(status).to_string(index=False), flush=True)

if __name__ == '__main__':
    main()
