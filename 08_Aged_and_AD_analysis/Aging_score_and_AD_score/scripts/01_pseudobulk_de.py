from pathlib import Path
import json
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from statsmodels.stats.multitest import multipletests
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# Parameters
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data/h5ad'
OUT = ROOT / 'outputs/01_region_pseudobulk'
GROUP_MAP = {'G21-4': 'Young', 'G21-7': 'Young', 'G21-1': 'Aged', 'G21-8': 'Aged', 'G21-2': 'AD', 'G21-5': 'AD'}
REGIONS = ['CTX_1', 'CTX_2', 'CTX_3', 'HIP', 'CA1', 'CA3', 'DG', 'TH', 'HY', 'STR', 'FT', 'VS']
MIN_COUNT, MIN_SAMPLES, FDR, MIN_LFC = 10, 2, 0.05, 0.25

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    # Shared genes
    features = []
    for sample in GROUP_MAP:
        a = ad.read_h5ad(DATA / f'{sample}_regions.h5ad', backed='r')
        assert a.var_names.is_unique
        features.append(set(a.var_names))
        if sample == 'G21-4':
            var = a.var[['Symbol', 'Gene_type']].copy()
        a.file.close()
    common = sorted(set.intersection(*features))
    var = var.loc[common].rename_axis('gene_id')
    var.to_csv(OUT / 'common_genes.csv')
    union = sorted(set.union(*features))
    pd.DataFrame({s: [g in f for g in union] for s, f in zip(GROUP_MAP, features)}, index=pd.Index(union, name='gene_id')).to_csv(OUT / 'gene_presence.csv.gz')
    blocks, meta, qc, means, obs_blocks = [], [], [], [], []
    for sample, condition in GROUP_MAP.items():
        a = ad.read_h5ad(DATA / f'{sample}_regions.h5ad')
        counts = sparse.csr_matrix(a.layers['counts'])
        assert np.isfinite(counts.data).all() and (counts.data >= 0).all() and np.equal(counts.data, np.round(counts.data)).all()
        x = sparse.csr_matrix(a.X)
        assert np.isfinite(x.data).all() and (x.data >= 0).all()
        idx = a.var_names.get_indexer(common)
        c = counts[:, idx]
        norm = x[:, idx]
        means.append(np.asarray(norm.mean(axis=0)).ravel())
        o = a.obs[['Regions']].copy()
        o['sample'] = sample
        o['condition'] = condition
        o['spot_id'] = a.obs_names.astype(str)
        o['x'], o['y'] = a.obsm['spatial'][:, 0], a.obsm['spatial'][:, 1]
        o['total_counts'] = np.asarray(counts.sum(axis=1)).ravel()
        o['n_genes'] = np.asarray((counts > 0).sum(axis=1)).ravel()
        obs_blocks.append(o.reset_index(drop=True))
        for region in sorted(o.Regions.unique()):
            mask = o.Regions.eq(region).to_numpy()
            row = {'sample_region': f'{sample}|{region}', 'sample': sample, 'condition': condition, 'region': region, 'n_spots': int(mask.sum()), 'in_DE_scope': region in REGIONS}
            qc.append({**row, 'total_counts': int(counts[mask].sum()), 'common_counts': int(c[mask].sum()), 'median_UMI': float(np.median(o.loc[mask, 'total_counts']))})
            if region in REGIONS:
                meta.append(row)
                blocks.append(np.asarray(c[mask].sum(axis=0)).ravel().astype(np.int64))
        print(f'Prepared {sample}: {a.shape}; common genes {len(common)}', flush=True)
        del a, counts, c, x, norm
    pb = pd.DataFrame(blocks, index=[r['sample_region'] for r in meta], columns=common)
    metadata = pd.DataFrame(meta).set_index('sample_region')
    pb.rename_axis('sample_region').to_csv(OUT / 'pseudobulk_counts.tsv.gz', sep='\t')
    metadata.to_csv(OUT / 'pseudobulk_metadata.csv')
    pd.DataFrame(qc).to_csv(OUT / 'sample_region_QC.csv', index=False)
    pd.DataFrame(means, index=pd.Index(GROUP_MAP, name='sample'), columns=common).to_csv(OUT / 'sample_mean_X.tsv.gz', sep='\t')
    pd.concat(obs_blocks, ignore_index=True).to_csv(OUT / 'spot_metadata.tsv.gz', sep='\t', index=False)
    # 同脑区六只小鼠共同拟合，提取两个对比
    results, model_qc = [], []
    for region in REGIONS:
        m = metadata.loc[metadata.region.eq(region)].copy()
        assert len(m) == 6 and (m.groupby('condition').size() == 2).all()
        c = pb.loc[m.index]
        keep = (c >= MIN_COUNT).sum(axis=0) >= MIN_SAMPLES
        c = c.loc[:, keep]
        print(f'Fitting {region}: {c.shape}', flush=True)
        dds = DeseqDataSet(counts=c, metadata=m, design='~condition', refit_cooks=True, n_cpus=4, quiet=True)
        dds.deseq2()
        n = pd.DataFrame(dds.layers['normed_counts'], index=m.index, columns=c.columns)
        n.rename_axis('sample_region').to_csv(OUT / f'{region}_normalized_pseudobulk.tsv.gz', sep='\t')
        model_qc.extend([{'region': region, 'sample': m.loc[k, 'sample'], 'size_factor': float(v), 'n_tested': c.shape[1]} for k, v in zip(m.index, dds.obs['size_factors'])])
        for condition in ['Aged', 'AD']:
            ds = DeseqStats(dds, contrast=['condition', condition, 'Young'], alpha=FDR, n_cpus=4, quiet=True)
            ds.summary()
            r = ds.results_df.copy().rename_axis('gene_id').reset_index()
            r['region'], r['comparison'] = region, f'{condition}_vs_Young'
            for s in GROUP_MAP:
                r[f'norm_{s}'] = n.loc[f'{s}|{region}'].reindex(r.gene_id).to_numpy()
            r['is_DEG'] = (r.padj < FDR) & (r.log2FoldChange.abs() > MIN_LFC)
            r['direction'] = np.where(r.log2FoldChange > 0, 'Up', 'Down')
            r.to_csv(OUT / f'{region}_{condition}_vs_Young.csv', index=False)
            results.append(r)
            print(f'{region} {condition}: {int(r.is_DEG.sum())} DEGs', flush=True)
    all_de = pd.concat(results, ignore_index=True).merge(var.reset_index(), on='gene_id', how='left')
    all_de['padj_across_regions'] = np.nan
    for comparison, ids in all_de.groupby('comparison').groups.items():
        ids = np.asarray(list(ids))
        valid = ids[all_de.loc[ids, 'pvalue'].notna().to_numpy()]
        all_de.loc[valid, 'padj_across_regions'] = multipletests(all_de.loc[valid, 'pvalue'], method='fdr_bh')[1]
    all_de.to_csv(OUT / 'all_region_DE_results.tsv.gz', sep='\t', index=False)
    all_de.loc[all_de.is_DEG].to_csv(OUT / 'DEGs_filtered.csv', index=False)
    all_de.groupby(['comparison', 'region', 'direction']).is_DEG.sum().rename('n_DEG').reset_index().to_csv(OUT / 'DEG_counts.csv', index=False)
    pd.DataFrame(model_qc).to_csv(OUT / 'model_QC.csv', index=False)
    (OUT / 'parameters.json').write_text(json.dumps({'groups': GROUP_MAP, 'regions': REGIONS, 'min_count': MIN_COUNT, 'min_samples': MIN_SAMPLES, 'fdr_within_region': FDR, 'min_abs_log2FC': MIN_LFC, 'n_common_genes': len(common), 'n_spots': sum(len(o) for o in obs_blocks), 'design': '~condition; six mice per region'}, indent=2))
    print('Pseudobulk DE completed', flush=True)

if __name__ == '__main__':
    main()
