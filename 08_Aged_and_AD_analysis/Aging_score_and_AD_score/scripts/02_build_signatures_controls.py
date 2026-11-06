from pathlib import Path
import json
import numpy as np
import pandas as pd

# Parameters
ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'data/references'
PB = ROOT / 'outputs/01_region_pseudobulk'
SIG = OUT = ROOT / 'outputs/02_signatures'
GROUP_MAP = {'G21-4': 'Young', 'G21-7': 'Young', 'G21-1': 'Aged', 'G21-8': 'Aged', 'G21-2': 'AD', 'G21-5': 'AD'}
REGIONS = ['CTX_1', 'CTX_2', 'CTX_3', 'HIP', 'CA1', 'CA3', 'DG', 'TH', 'HY', 'STR', 'FT', 'VS']
NBIN, CTRL, SEED = 24, 100, 1

def direction_audit(rows):
    g = rows.groupby('gene_id')
    result = g.agg(symbol=('symbol', 'first'), n_records=('direction', 'size'), n_directions=('direction', 'nunique'), direction=('direction', 'first')).reset_index()
    result['conflict'] = result.n_directions > 1
    result.loc[result.conflict, 'direction'] = 'Conflict'
    return result

def main():
    for d in [SIG, OUT]:
        d.mkdir(parents=True, exist_ok=True)
    var = pd.read_csv(PB / 'common_genes.csv').set_index('gene_id')
    genes = var.index.to_numpy()
    # 外部方向，不做层级计票
    allen = pd.read_csv(BASE / 'Allen_aging_degenes.csv')
    allen = allen.loc[allen.grouping_type.isin(['subclass', 'supertype']) & (allen.adjusted_pvalue < 0.01)].copy()
    allen = allen.rename(columns={'gene_identifier': 'gene_id', 'gene_symbol': 'symbol'})
    allen['direction'] = np.where(allen.age_effect_size > 0, 'Up', 'Down')
    allen.to_csv(SIG / 'Aging_external_records.csv', index=False)
    for level, frame in allen.groupby('grouping_type'):
        direction_audit(frame).to_csv(SIG / f'Aging_external_{level}_audit.csv', index=False)
    ad_rows = pd.read_excel(BASE / '2023-NG-AD.xlsx', sheet_name='Supplementary Table 12', header=1)
    ad_rows = ad_rows.loc[ad_rows['Q-value'] < 0.05].copy()
    ad_rows['gene_id'] = ad_rows['Gene ID'].str.replace(r'\.\d+$', '', regex=True)
    ad_rows['symbol'] = ad_rows['Gene name']
    assert set(ad_rows['Max condition']) == {'Wild type', 'EOAD (5xFAD)'}
    ad_rows['direction'] = np.where(ad_rows['Max condition'].eq('EOAD (5xFAD)'), 'Up', 'Down')
    ad_rows.to_csv(SIG / 'AD_external_records.csv', index=False)
    de = pd.read_csv(PB / 'all_region_DE_results.tsv.gz', sep='\t')
    signatures, flow = {}, []
    for score, rows, contrast in [('Aging', allen, 'Aged_vs_Young'), ('AD', ad_rows, 'AD_vs_Young')]:
        ext = direction_audit(rows)
        ext['in_common_genes'] = ext.gene_id.isin(genes)
        ext.to_csv(SIG / f'{score}_external_direction_audit.csv', index=False)
        d = de.loc[de.comparison.eq(contrast)].copy()
        hits = d.loc[d.is_DEG].rename(columns={'Symbol': 'symbol'})
        own = direction_audit(hits)
        supports = hits.groupby('gene_id').agg(support_regions=('region', lambda s: ';'.join(sorted(s))), min_padj=('padj', 'min'), min_padj_across_regions=('padj_across_regions', 'min'))
        own = own.merge(supports, on='gene_id', how='left')
        own.to_csv(SIG / f'{score}_own_direction_audit.csv', index=False)
        audit = ext.merge(own.drop(columns='symbol'), on='gene_id', how='outer', suffixes=('_external', '_own'))
        audit['selected'] = audit.direction_external.isin(['Up', 'Down']) & audit.direction_external.eq(audit.direction_own) & audit.gene_id.isin(genes)
        audit['reason'] = np.select([audit.selected, audit.direction_external.eq('Conflict'), audit.direction_own.eq('Conflict'), audit.direction_external.isna(), audit.direction_own.isna()], ['selected', 'external_direction_conflict', 'own_region_direction_conflict', 'no_external_evidence', 'no_own_DEG'], default='direction_mismatch')
        audit.to_csv(SIG / f'{score}_intersection_audit.csv', index=False)
        sig = audit.loc[audit.selected, ['gene_id', 'symbol', 'direction_external', 'support_regions', 'min_padj', 'min_padj_across_regions']].rename(columns={'direction_external': 'direction'}).sort_values(['direction', 'gene_id'])
        sig['passes_global_sensitivity'] = sig.min_padj_across_regions < 0.05
        sig.to_csv(SIG / f'{score}_signature.csv', index=False)
        signatures[score] = sig
        flow.append({'score': score, 'external_genes': len(ext), 'external_conflict': int(ext.conflict.sum()), 'external_nonconflict_common': int((~ext.conflict & ext.in_common_genes).sum()), 'own_DEG_unique': len(own), 'own_conflict': int(own.conflict.sum()), 'final_Up': int(sig.direction.eq('Up').sum()), 'final_Down': int(sig.direction.eq('Down').sum()), 'global_sensitivity_retained': int(sig.passes_global_sensitivity.sum())})
        print(flow[-1], flush=True)
    pd.DataFrame(flow).to_csv(SIG / 'signature_flow.csv', index=False)
    # 固定背景：六只小鼠均值等权分箱，各样本共用 control genes
    means = pd.read_csv(PB / 'sample_mean_X.tsv.gz', sep='\t', index_col=0).reindex(columns=genes).mean(axis=0).to_numpy()
    rng = np.random.default_rng(SEED)
    order = np.argsort(means + rng.normal(size=len(means)) * 1e-30)
    bins = np.empty(len(genes), int)
    bins[order] = np.minimum(np.arange(len(genes)) * NBIN // len(genes), NBIN - 1)
    pd.DataFrame({'gene_id': genes, 'sample_equal_mean_X': means, 'bin': bins}).to_csv(OUT / 'background_bins.csv', index=False)
    modules, control_rows, status = {}, [], []
    for name, sig in signatures.items():
        sizes = sig.direction.value_counts()
        available = sizes.get('Up', 0) > 0 and sizes.get('Down', 0) > 0
        status.append({'score': name, 'available': available, 'Up': int(sizes.get('Up', 0)), 'Down': int(sizes.get('Down', 0)), 'small_module': min(sizes.get('Up', 0), sizes.get('Down', 0)) < 5})
        if not available:
            continue
        for direction in ['Up', 'Down']:
            idx = var.index.get_indexer(sig.loc[sig.direction.eq(direction), 'gene_id'])
            controls = set()
            for i in idx:
                pool = np.flatnonzero(bins == bins[i])
                controls.update(rng.choice(pool, min(CTRL, len(pool)), replace=False).tolist())
            ctrl = np.array(sorted(controls))
            modules[f'{name}_{direction}'] = (idx, ctrl)
            control_rows.extend({'module': f'{name}_{direction}', 'gene_id': genes[i]} for i in ctrl)
    pd.DataFrame(status).to_csv(OUT / 'score_status.csv', index=False)
    pd.DataFrame(control_rows, columns=['module', 'gene_id']).to_csv(OUT / 'module_controls.csv', index=False)
    var.reset_index().to_csv(OUT / 'common_genes.csv', index=False)
    print('Signatures and fixed controls saved:', OUT)

if __name__ == '__main__':
    main()
