#!/usr/bin/env python3
"""
step0_inspect_adatas_regions.py
================================
检查 6 个样本的新 h5ad 文件，为后续 rMATS 脚本提供依据。

当前版本适配：
  - h5ad 目录：/path/to/project/20260705_Fig3/anno_h5ads
  - h5ad 文件名：<SAMPLE>_regions.h5ad
  - 区域注释列：adata.obs["Regions"]
  - rMATS 工作目录：/path/to/project/20260705_Fig3/rmats_works

用法：
  python step0_inspect_adatas_regions.py

检查项目：
  1. 文件是否存在 & 基本信息（n_obs, n_vars）
  2. obs_names 格式（前 5 个）及 barcode 长度分布
  3. obs 中关键列是否存在：Regions、Spot_col、Spot_row、sample_id
  4. Regions 分布：包含排除区域与有效区域
  5. sample_id 分布
  6. 样本 × Regions spot 数量矩阵
  7. 跨样本 Regions 名称一致性检查

输出：
  <WORK_DIR>/metadata/step0_regions_spot_counts_all.tsv
  <WORK_DIR>/metadata/step0_regions_spot_counts_valid.tsv
  <WORK_DIR>/metadata/step0_regions_presence_valid.tsv
  <WORK_DIR>/metadata/step0_sample_summary.tsv
"""

from __future__ import annotations

import collections
import os
import sys
from typing import Dict, Optional, Set

import pandas as pd

# ============================================================
# Configuration
# ============================================================
ADATA_DIR = "/path/to/project/20260705_Fig3/anno_h5ads"
WORK_DIR = "/path/to/project/20260705_Fig3/rmats_works"
META_DIR = os.path.join(WORK_DIR, "metadata")

REGION_COL = "Regions"
EXCLUDE_CLUSTERS = {"unknown", "EPI", "CTX_4"}

GROUPS = {
    "control": ["G21-4", "G21-7"],
    "old": ["G21-1", "G21-8"],
    "AD": ["G21-2", "G21-5"],
}

SAMPLES = ["G21-1", "G21-2", "G21-4", "G21-5", "G21-7", "G21-8"]
H5AD_SUFFIX = "_regions.h5ad"

REQUIRED_COLS = [REGION_COL, "Spot_col", "Spot_row"]
OPTIONAL_COLS = ["sample_id", "reg", "brain_cluster_fullname", "brain_cluster_shortname"]

# ============================================================
# 导入 scanpy
# ============================================================
try:
    import scanpy as sc
except ImportError:
    print("[ERROR] 请先安装 scanpy：pip install scanpy")
    sys.exit(1)


def pure_barcode(obs_name: str, sample: str) -> str:
    """从 obs_name 中去掉 _<SAMPLE> 后缀，得到 Space Ranger 原始 barcode。"""
    sample_tag = f"_{sample}"
    if obs_name.endswith(sample_tag):
        return obs_name[: -len(sample_tag)]
    return obs_name


def safe_value_set(series: pd.Series) -> Set[str]:
    """将区域列转为字符串集合，同时去掉缺失值。"""
    return set(series.dropna().astype(str).tolist())


print("=" * 75)
print("Step0: h5ad 文件检查（Regions 版本）")
print("=" * 75)
print(f"h5ad 数据目录: {ADATA_DIR}")
print(f"rMATS 工作目录: {WORK_DIR}")
print(f"区域列: {REGION_COL}")
print(f"排除区域: {sorted(EXCLUDE_CLUSTERS)}")
print("样本分组:")
for grp, slist in GROUPS.items():
    print(f"  {grp:10s}: {slist}")
print()

os.makedirs(META_DIR, exist_ok=True)

all_regions: Dict[str, Set[str]] = {}
valid_regions_by_sample: Dict[str, Set[str]] = {}
summary_rows = []
count_rows_all = []
count_rows_valid = []

# ============================================================
# 逐样本检查
# ============================================================
for sample in SAMPLES:
    h5ad = os.path.join(ADATA_DIR, f"{sample}{H5AD_SUFFIX}")

    print("=" * 75)
    print(f"样本: {sample}")
    print(f"  路径: {h5ad}")

    if not os.path.exists(h5ad):
        print("  [ERROR] 文件不存在！")
        summary_rows.append({
            "sample": sample,
            "h5ad_path": h5ad,
            "status": "missing_file",
            "n_obs": None,
            "n_vars": None,
            "n_regions_all": None,
            "n_regions_valid": None,
        })
        all_regions[sample] = set()
        valid_regions_by_sample[sample] = set()
        continue

    try:
        adata = sc.read(h5ad)
    except Exception as e:
        print(f"  [ERROR] 读取失败: {e}")
        summary_rows.append({
            "sample": sample,
            "h5ad_path": h5ad,
            "status": f"read_error: {e}",
            "n_obs": None,
            "n_vars": None,
            "n_regions_all": None,
            "n_regions_valid": None,
        })
        all_regions[sample] = set()
        valid_regions_by_sample[sample] = set()
        continue

    print(f"  n_obs={adata.n_obs}, n_vars={adata.n_vars}")

    # --------------------------------------------------------
    # obs_names / barcode 格式
    # --------------------------------------------------------
    print("\n  [1] obs_names 格式（前 5 个）:")
    for name in adata.obs_names[:5]:
        print(f"      '{name}'")

    print("\n  [2] barcode 长度分析:")
    lengths = [len(pure_barcode(str(name), sample)) for name in adata.obs_names]
    len_dist = collections.Counter(lengths)
    for length, cnt in sorted(len_dist.items()):
        flag = "✓" if length == 16 else "⚠"
        print(f"      {flag} 长度 {length}: {cnt} 个 barcode")

    # --------------------------------------------------------
    # obs 关键列检查
    # --------------------------------------------------------
    print("\n  [3] obs 关键列检查:")
    missing_required = []
    for col in REQUIRED_COLS:
        if col in adata.obs.columns:
            print(f"      ✓  {col}")
        else:
            print(f"      ✗ 缺失  {col}")
            missing_required.append(col)

    print("      --- 可选列 ---")
    for col in OPTIONAL_COLS:
        flag = "✓" if col in adata.obs.columns else "-"
        print(f"      {flag}  {col}")

    print("      --- obs 全部列 ---")
    for col in adata.obs.columns:
        dtype = str(adata.obs[col].dtype)
        n_unique = adata.obs[col].nunique(dropna=True)
        print(f"           {col:<30} dtype={dtype:<15} nunique={n_unique}")

    if missing_required:
        print(f"\n  [ERROR] 缺失必要列 {missing_required}，跳过该样本的 Regions 统计。")
        summary_rows.append({
            "sample": sample,
            "h5ad_path": h5ad,
            "status": "missing_required_columns",
            "n_obs": adata.n_obs,
            "n_vars": adata.n_vars,
            "n_regions_all": None,
            "n_regions_valid": None,
        })
        all_regions[sample] = set()
        valid_regions_by_sample[sample] = set()
        print()
        continue

    # --------------------------------------------------------
    # Regions 分布
    # --------------------------------------------------------
    region_series = adata.obs[REGION_COL].astype(str)
    vc = region_series.value_counts(dropna=False)
    regions_this = set(vc.index.astype(str))
    valid_regions_this = regions_this - EXCLUDE_CLUSTERS

    all_regions[sample] = regions_this
    valid_regions_by_sample[sample] = valid_regions_this

    print(f"\n  [4] {REGION_COL} 分布:")
    for region, cnt in vc.items():
        region_str = str(region)
        flag = "✗(排除)" if region_str in EXCLUDE_CLUSTERS else "✓"
        print(f"      {flag:<8} {region_str:<20} {int(cnt):>5} spots")

    print(f"\n  有效区域数: {len(valid_regions_this)}")
    print(f"  有效区域: {sorted(valid_regions_this)}")

    # --------------------------------------------------------
    # sample_id 分布
    # --------------------------------------------------------
    if "sample_id" in adata.obs.columns:
        print("\n  [5] sample_id 分布:")
        print(adata.obs["sample_id"].value_counts(dropna=False).to_string())
    else:
        print("\n  [5] sample_id 列: 不存在")

    # --------------------------------------------------------
    # 汇总行
    # --------------------------------------------------------
    row_all = {"sample": sample}
    row_all.update({str(k): int(v) for k, v in vc.to_dict().items()})
    count_rows_all.append(row_all)

    row_valid = {"sample": sample}
    for region, cnt in vc.to_dict().items():
        region_str = str(region)
        if region_str not in EXCLUDE_CLUSTERS:
            row_valid[region_str] = int(cnt)
    count_rows_valid.append(row_valid)

    summary_rows.append({
        "sample": sample,
        "h5ad_path": h5ad,
        "status": "ok",
        "n_obs": adata.n_obs,
        "n_vars": adata.n_vars,
        "n_regions_all": len(regions_this),
        "n_regions_valid": len(valid_regions_this),
        "regions_all": ";".join(sorted(regions_this)),
        "regions_valid": ";".join(sorted(valid_regions_this)),
    })

    print()

# ============================================================
# 跨样本 Regions 一致性检查
# ============================================================
print("=" * 75)
print(f"跨样本 {REGION_COL} 一致性检查（已排除 {sorted(EXCLUDE_CLUSTERS)}）")
print("=" * 75)

valid_samples = [s for s in SAMPLES if valid_regions_by_sample.get(s)]
if len(valid_samples) < 2:
    print("  有效样本数不足，跳过一致性检查")
else:
    union_regions = set.union(*[valid_regions_by_sample[s] for s in valid_samples])
    inter_regions = set.intersection(*[valid_regions_by_sample[s] for s in valid_samples])

    print(f"  所有有效区域并集：{len(union_regions)} 个")
    print(f"  所有样本共有有效区域：{len(inter_regions)} 个")
    for region in sorted(inter_regions):
        print(f"    ✓  {region}")

    only_in_some = union_regions - inter_regions
    if only_in_some:
        print(f"\n  仅部分样本有的有效区域：{len(only_in_some)} 个")
        for region in sorted(only_in_some):
            present_in = [s for s in valid_samples if region in valid_regions_by_sample.get(s, set())]
            missing_in = [s for s in valid_samples if region not in valid_regions_by_sample.get(s, set())]
            print(f"    ⚠  {region:<20} 出现于: {present_in}; 缺失于: {missing_in}")
    else:
        print("\n  ✓ 所有样本有效区域完全一致，无缺失")

# ============================================================
# 输出样本 × Regions spot 数量矩阵
# ============================================================
print()
print("=" * 75)
print("样本 × Regions spot 数量汇总")
print("=" * 75)

summary_df = pd.DataFrame(summary_rows)
summary_out = os.path.join(META_DIR, "step0_sample_summary.tsv")
summary_df.to_csv(summary_out, sep="\t", index=False)

if count_rows_all:
    df_all = pd.DataFrame(count_rows_all).set_index("sample").fillna(0).astype(int)
    cols = sorted([c for c in df_all.columns if c not in EXCLUDE_CLUSTERS])
    excluded_cols = sorted([c for c in df_all.columns if c in EXCLUDE_CLUSTERS])
    df_all = df_all[cols + excluded_cols]
    all_out = os.path.join(META_DIR, "step0_regions_spot_counts_all.tsv")
    df_all.to_csv(all_out, sep="\t")
    print("\n[全部区域；排除区域放在最后]")
    print(df_all.to_string())
    print(f"\n  → {all_out}")

if count_rows_valid:
    df_valid = pd.DataFrame(count_rows_valid).set_index("sample").fillna(0).astype(int)
    df_valid = df_valid[sorted(df_valid.columns)]
    valid_out = os.path.join(META_DIR, "step0_regions_spot_counts_valid.tsv")
    df_valid.to_csv(valid_out, sep="\t")
    print("\n[有效区域；已排除 unknown/EPI/CTX_4]")
    print(df_valid.to_string())
    print(f"\n  → {valid_out}")

    presence_df = (df_valid > 0).astype(int)
    presence_out = os.path.join(META_DIR, "step0_regions_presence_valid.tsv")
    presence_df.to_csv(presence_out, sep="\t")
    print(f"  → {presence_out}")

print(f"  → {summary_out}")

print()
print("=" * 75)
print("检查完成。")
print("=" * 75)
