#!/usr/bin/env python3
"""
step1_extract_barcodes_regions.py
=================================
从新 h5ad 文件中按 adata.obs["Regions"] 提取各脑区 barcode 列表。

当前版本适配：
  - h5ad 文件名：<SAMPLE>_regions.h5ad
  - 区域注释列：adata.obs["Regions"]
  - rMATS 工作目录：/path/to/project/20260705_Fig3/rmats_works
  - 排除区域：unknown、EPI、CTX_4

用法：
  python step1_extract_barcodes_regions.py <样本名> <h5ad路径>

示例：
  python step1_extract_barcodes_regions.py G21-4 \
      /path/to/project/20260705_Fig3/anno_h5ads/G21-4_regions.h5ad

输出目录：
  <WORK_DIR>/barcodes/<SAMPLE>/
  ├── <SAMPLE>_barcode_to_pixel.tsv      # barcode + Spot_col + Spot_row + Regions
  ├── <SAMPLE>_barcode_to_region.tsv     # barcode → Regions（无表头）
  ├── <SAMPLE>_region_list.tsv           # 有效区域列表：region_safe_name\tregion_original_name
  ├── <SAMPLE>_ALL_barcodes.txt          # 全脑 barcode（排除 unknown/EPI/CTX_4）
  └── <SAMPLE>_<REGION>_barcodes.txt     # 各区域 barcode（samtools -D CB: 用）

说明：
  当前 Regions 名称均为无空格短名（如 CTX_1、STR、HIP），region_safe_name 与原名通常一致。
  为了后续脚本稳健性，仍保留 safe_name 机制，避免未来区域名出现空格或特殊字符时破坏文件名。
"""

from __future__ import annotations

import os
import re
import sys
from typing import Dict, Iterable

import pandas as pd
import scanpy as sc

# ============================================================
# Configuration
# ============================================================
WORK_DIR = "/path/to/project/20260705_Fig3/rmats_works"
BC_BASE = os.path.join(WORK_DIR, "barcodes")
LOG_DIR = os.path.join(WORK_DIR, "logs")

REGION_COL = "Regions"
EXCLUDE_CLUSTERS = {"unknown", "EPI", "CTX_4"}

REQUIRED_COLS = [REGION_COL, "Spot_col", "Spot_row"]


def make_safe_name(value: object) -> str:
    """生成适合放入文件名和 bash 数组的 region 名称。"""
    text = str(value).strip()
    text = re.sub(r"\s+", "_", text)
    text = re.sub(r"[^0-9A-Za-z._-]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text if text else "NA"


def pure_barcode(obs_name: str, sample: str) -> str:
    """从 obs_name 中去掉 _<SAMPLE> 后缀，得到 Space Ranger 原始 barcode。"""
    sample_tag = f"_{sample}"
    if obs_name.endswith(sample_tag):
        return obs_name[: -len(sample_tag)]
    return obs_name


def write_lines(path: str, values: Iterable[str]) -> None:
    values = list(values)
    with open(path, "w") as f:
        if values:
            f.write("\n".join(values) + "\n")
        else:
            f.write("")


# ============================================================
# Parse arguments
# ============================================================
if len(sys.argv) != 3:
    print(f"用法: python {sys.argv[0]} <样本名> <h5ad路径>")
    sys.exit(1)

SAMPLE = sys.argv[1]
H5AD_PATH = sys.argv[2]
BC_DIR = os.path.join(BC_BASE, SAMPLE)

# ============================================================
# Check inputs and create output directories
# ============================================================
if not os.path.exists(H5AD_PATH):
    print(f"[ERROR] h5ad 文件不存在: {H5AD_PATH}")
    sys.exit(1)

os.makedirs(BC_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)

print(f"[Step1] 样本: {SAMPLE}")
print(f"  h5ad:       {H5AD_PATH}")
print(f"  工作目录:   {WORK_DIR}")
print(f"  输出目录:   {BC_DIR}")
print(f"  区域列:     {REGION_COL}")
print(f"  排除区域:   {sorted(EXCLUDE_CLUSTERS)}")
print()

# ============================================================
# Read AnnData
# ============================================================
try:
    adata = sc.read(H5AD_PATH)
except Exception as e:
    print(f"[ERROR] 读取 h5ad 失败: {e}")
    sys.exit(1)

print(f"  总 spot 数: {adata.n_obs}")

# ============================================================
# 检查必要列
# ============================================================
missing_cols = [col for col in REQUIRED_COLS if col not in adata.obs.columns]
if missing_cols:
    print(f"[ERROR] adata.obs 缺少必要列: {missing_cols}")
    print("        当前 obs.columns:")
    for col in adata.obs.columns:
        print(f"          - {col}")
    sys.exit(1)

# ============================================================
# 提取纯 barcode（去掉 _<SAMPLE> 后缀）
# ============================================================
adata.obs = adata.obs.copy()
adata.obs["barcode"] = [pure_barcode(str(idx), SAMPLE) for idx in adata.obs_names]

sample_bc = str(adata.obs["barcode"].iloc[0])
print(f"  barcode 示例: '{sample_bc}'（长度 {len(sample_bc)}）")
if len(sample_bc) != 16:
    print("  ⚠  barcode 长度不是 16 bp，请检查 obs_names 格式！")

# ============================================================
# 标准化 Regions 列为字符串
# ============================================================
adata.obs[REGION_COL] = adata.obs[REGION_COL].astype(str)
all_regions = adata.obs[REGION_COL].dropna().unique().tolist()

# safe name 映射；如果出现 safe name 冲突，直接报错，避免覆盖 barcode 文件
safe_to_original: Dict[str, str] = {}
original_to_safe: Dict[str, str] = {}
for region in sorted(all_regions):
    safe = make_safe_name(region)
    if safe in safe_to_original and safe_to_original[safe] != region:
        print("[ERROR] 不同 Regions 名称生成了相同 safe name，可能导致文件覆盖：")
        print(f"        {safe_to_original[safe]}  vs  {region}  →  {safe}")
        sys.exit(1)
    safe_to_original[safe] = region
    original_to_safe[region] = safe

valid_regions = sorted([region for region in all_regions if region not in EXCLUDE_CLUSTERS])

print("\n  脑区统计:")
vc = adata.obs[REGION_COL].value_counts(dropna=False)
for region in sorted(all_regions):
    flag = "✓" if region not in EXCLUDE_CLUSTERS else "✗(排除)"
    safe = original_to_safe[region]
    print(f"    {flag:<8} {region:<20} {int(vc.get(region, 0)):>5} spots  safe_name={safe}")
print(f"\n  有效脑区（{len(valid_regions)} 个）: {valid_regions}")
print()

# ============================================================
# 输出1: barcode_to_pixel.tsv（空间可视化用）
# ============================================================
mask_valid = ~adata.obs[REGION_COL].isin(EXCLUDE_CLUSTERS)

pixel_df = adata.obs.loc[mask_valid, ["barcode", "Spot_col", "Spot_row", REGION_COL]].copy()
pixel_out = os.path.join(BC_DIR, f"{SAMPLE}_barcode_to_pixel.tsv")
pixel_df.to_csv(pixel_out, sep="\t", index=False)
print(f"  → {os.path.basename(pixel_out):<35} ({len(pixel_df)} spots)")

# ============================================================
# 输出2: barcode_to_region.tsv（无表头，barcode\tRegions）
# ============================================================
region_df = adata.obs.loc[mask_valid, ["barcode", REGION_COL]].copy()
region_out = os.path.join(BC_DIR, f"{SAMPLE}_barcode_to_region.tsv")
region_df.to_csv(region_out, sep="\t", index=False, header=False)
print(f"  → {os.path.basename(region_out):<35} ({len(region_df)} spots)")

# 为兼容后续旧脚本，也额外输出旧命名 barcode_to_cluster.tsv
cluster_compat_out = os.path.join(BC_DIR, f"{SAMPLE}_barcode_to_cluster.tsv")
region_df.to_csv(cluster_compat_out, sep="\t", index=False, header=False)
print(f"  → {os.path.basename(cluster_compat_out):<35} ({len(region_df)} spots; compat)")

# ============================================================
# 输出3: region_list.tsv（有效脑区列表）
# 第一列 safe name，用于文件名和后续 bash；第二列原始 Regions 名称
# ============================================================
region_list_out = os.path.join(BC_DIR, f"{SAMPLE}_region_list.tsv")
with open(region_list_out, "w") as f:
    for region in valid_regions:
        f.write(f"{original_to_safe[region]}\t{region}\n")
print(f"  → {os.path.basename(region_list_out):<35} ({len(valid_regions)} 个有效脑区)")

# 为兼容当前 Step2 旧脚本，也额外输出旧命名 cluster_list.tsv
cluster_list_compat_out = os.path.join(BC_DIR, f"{SAMPLE}_cluster_list.tsv")
with open(cluster_list_compat_out, "w") as f:
    for region in valid_regions:
        f.write(f"{original_to_safe[region]}\t{region}\n")
print(f"  → {os.path.basename(cluster_list_compat_out):<35} ({len(valid_regions)} 个有效脑区; compat)")

# ============================================================
# 输出4: ALL barcodes（全脑，排除 EXCLUDE_CLUSTERS）
# ============================================================
barcodes_all = adata.obs.loc[mask_valid, "barcode"].astype(str).tolist()
out_all = os.path.join(BC_DIR, f"{SAMPLE}_ALL_barcodes.txt")
write_lines(out_all, barcodes_all)
print(f"  → {os.path.basename(out_all):<35} ({len(barcodes_all)} spots)")

# ============================================================
# 输出5: 各脑区 barcode 文件（samtools -D CB: 过滤用）
# 文件名使用 safe name；筛选仍使用原始 Regions 名称
# ============================================================
print("\n  各脑区 barcode 文件:")
for region in valid_regions:
    safe = original_to_safe[region]
    mask = adata.obs[REGION_COL] == region
    barcodes = adata.obs.loc[mask, "barcode"].astype(str).tolist()
    out_file = os.path.join(BC_DIR, f"{SAMPLE}_{safe}_barcodes.txt")
    write_lines(out_file, barcodes)
    print(f"    {safe:<20} {len(barcodes):>5} spots → {os.path.basename(out_file)}")

# ============================================================
# 输出6: 运行摘要
# ============================================================
summary_df = pd.DataFrame({
    "sample": [SAMPLE],
    "h5ad_path": [H5AD_PATH],
    "work_dir": [WORK_DIR],
    "barcode_dir": [BC_DIR],
    "region_col": [REGION_COL],
    "exclude_regions": [";".join(sorted(EXCLUDE_CLUSTERS))],
    "n_obs_total": [adata.n_obs],
    "n_spots_valid": [int(mask_valid.sum())],
    "n_regions_all": [len(all_regions)],
    "n_regions_valid": [len(valid_regions)],
    "regions_valid": [";".join(valid_regions)],
})
summary_out = os.path.join(BC_DIR, f"{SAMPLE}_step1_summary.tsv")
summary_df.to_csv(summary_out, sep="\t", index=False)
print(f"\n  → {os.path.basename(summary_out):<35} (summary)")

print()
print(f"[Step1] 完成 → {BC_DIR}")
print(f"  下一步示例: bash step2_generate_pseudobulk_bam_regions.sh {SAMPLE} <BAM路径>")
