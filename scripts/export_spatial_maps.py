"""
从 h5ad 导出 RNA 编辑流程所需的空间映射文件
输出（每个样本各一套）：
  cb_to_pixel_{sample}.tsv         格式: AACGCCTCAATCTCGA\t139x61
  pixel_to_clu_ID_{sample}.tsv     格式: 139x61\t2
  clu_ID_to_region_{sample}.tsv    格式: 2\tCortex_1
"""

import anndata as ad
import pandas as pd
import os

SAMPLES = {
    # 'G21-4':       '/histor/zhao/pangkun/Zhaolab/Lupeng/Zhujunjie/20260427_AS_analysis/adatas/G21-4_anno.h5ad',
    # 'G21-7':       '/histor/zhao/pangkun/Zhaolab/Lupeng/Zhujunjie/20260427_AS_analysis/adatas/G21-7_anno.h5ad',
    'ZZJ_Reg1_10': '/histor/zhao/pangkun/Zhaolab/Lupeng/Zhujunjie/20260307_testdata/pangkun_st_data/ZZJ_Reg1_10_anno.h5ad',
    'WT1_10':      '/histor/zhao/pangkun/Zhaolab/Lupeng/Zhujunjie/20260307_testdata/pangkun_st_data/WT1_10_anno.h5ad',
}

# 各样本的坐标列名（col, row）
COL_MAP = {
    'G21-4':       ('Spot_col', 'Spot_row'),
    'G21-7':       ('Spot_col', 'Spot_row'),
    'ZZJ_Reg1_10': ('Spot_col', 'Spot_row'),
    'WT1_10':      ('array_col', 'array_row'),  # 10x Visium
}

OUT_DIR = '/histor/zhao/pangkun/Zhaolab/Lupeng/Zhujunjie/20260505_AtoI/data'
os.makedirs(OUT_DIR, exist_ok=True)

for sample_id, h5ad_path in SAMPLES.items():
    print(f'\n=== {sample_id} ===')
    adata = ad.read_h5ad(h5ad_path)
    obs = adata.obs.copy()

    col_col, row_col = COL_MAP[sample_id]

    # CB: obs_names 格式 AACGCCTCAATCTCGA_G21-4 或裸CB，统一取下划线前部分
    obs['CB'] = obs.index.str.split('_').str[0]

    # pixel: colxrow
    obs['pixel'] = obs[col_col].astype(str) + 'x' + obs[row_col].astype(str)

    obs['region'] = obs['AS_cluster']
    obs['clu']    = obs['leiden'].astype(str)

    # ── 文件1: cb_to_pixel ──────────────────────────────────────
    cb_to_pixel = obs[['CB', 'pixel']].drop_duplicates()
    out1 = os.path.join(OUT_DIR, f'cb_to_pixel_{sample_id}.tsv')
    cb_to_pixel.to_csv(out1, sep='\t', index=False, header=False)
    print(f'  cb_to_pixel       -> {out1}  ({len(cb_to_pixel)} barcodes)')

    # ── 文件2: pixel_to_clu_ID ──────────────────────────────────
    pixel_to_clu = obs[['pixel', 'clu']].drop_duplicates()
    out2 = os.path.join(OUT_DIR, f'pixel_to_clu_ID_{sample_id}.tsv')
    pixel_to_clu.to_csv(out2, sep='\t', index=False, header=False)
    print(f'  pixel_to_clu_ID   -> {out2}  ({len(pixel_to_clu)} pixels)')

    # ── 文件3: clu_ID_to_region ─────────────────────────────────
    clu_to_region = obs[['clu', 'region']].drop_duplicates().sort_values('clu')
    out3 = os.path.join(OUT_DIR, f'clu_ID_to_region_{sample_id}.tsv')
    clu_to_region.to_csv(out3, sep='\t', index=False, header=False)
    print(f'  clu_ID_to_region  -> {out3}  ({len(clu_to_region)} clusters)')
    for _, row in clu_to_region.iterrows():
        print(f'    {row["clu"]:>4}  ->  {row["region"]}')

print('\n全部导出完成。')
