"""Prepare gene results and AS-testable events from DEG and rMATS JC tables.
Run scripts 00 -> 01 -> 02 -> 03. Requires numpy, pandas and matplotlib.
Inputs: data/spot_DEG_all_regions_all_comparisons.tsv.gz and
        data/rmats_results/{control_vs_AD,control_vs_old}/{region}/*.MATS.JC.txt.
rMATS sample group 1 is Young; group 2 is AD or Aged.
"""
from pathlib import Path
import numpy as np
import pandas as pd

# Parameters
DATA_DIR = Path(__file__).resolve().parent / "data"
OUT_DIR = Path(__file__).resolve().parent / "results"
GENE_FILE = DATA_DIR / "spot_DEG_all_regions_all_comparisons.tsv.gz"
RMATS_DIR = DATA_DIR / "rmats_results"
COMPARISONS = {"control_vs_AD": "AD_vs_Young", "control_vs_old": "Aged_vs_Young"}
REGIONS = ["CTX_1", "CTX_2", "CTX_3", "HIP", "CA1", "CA3", "DG", "TH", "HY", "STR", "FT", "VS"]
READ_THRESH = 10
EVENT_COLS = {
    "SE": ["chr", "strand", "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
    "MXE": ["chr", "strand", "1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
    "A3SS": ["chr", "strand", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"],
    "A5SS": ["chr", "strand", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"],
    "RI": ["chr", "strand", "riExonStart_0base", "riExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
}
OUT_DIR.mkdir(parents=True, exist_ok=True)

# 整理完整gene差异表，不能只输入显著DEG
# 必需列：comparison, region, Symbol, gene_id, log2FC, padj
gene = pd.read_csv(GENE_FILE, sep="\t")
gene = gene[gene.comparison.isin(COMPARISONS.values()) & gene.region.isin(REGIONS)].copy()
gene["gene_name"] = gene.Symbol.astype("string").str.strip().str.replace('"', '', regex=False)
gene = gene[gene.gene_name.notna() & gene.log2FC.notna()]
gene = gene.sort_values("padj", na_position="last").drop_duplicates(["region", "comparison", "gene_name"])
gene[["comparison", "region", "gene_name", "gene_id", "log2FC", "padj"]].to_csv(OUT_DIR / "gene_results.tsv.gz", sep="\t", index=False)

# 读取rMATS事件
def split_counts(x):
    return x.astype(str).str.split(",", expand=True).apply(pd.to_numeric, errors="coerce").fillna(0)

event_tables = []

for rmats_comparison, comparison in COMPARISONS.items():
    for region in REGIONS:
        for ast, coordinates in EVENT_COLS.items():
            file = RMATS_DIR / rmats_comparison / region / f"{ast}.MATS.JC.txt"
            usecols = list(dict.fromkeys([
                "GeneID", "geneSymbol", "IJC_SAMPLE_1", "SJC_SAMPLE_1",
                "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncLevelDifference", "FDR",
            ] + coordinates))
            df = pd.read_csv(file, sep="\t", usecols=usecols)

            total1 = split_counts(df.IJC_SAMPLE_1) + split_counts(df.SJC_SAMPLE_1)
            total2 = split_counts(df.IJC_SAMPLE_2) + split_counts(df.SJC_SAMPLE_2)
            coverage = (total1 >= READ_THRESH).all(axis=1) & (total2 >= READ_THRESH).all(axis=1)
            # Condition - Young；不预先筛选FDR或ΔPSI幅度
            dpsi = -pd.to_numeric(df.IncLevelDifference, errors="coerce")
            fdr = pd.to_numeric(df.FDR, errors="coerce")
            keep = coverage & dpsi.notna() & fdr.notna()
            tmp = df.loc[keep].copy()
            tmp["dPSI"] = dpsi.loc[keep]
            tmp["FDR_num"] = fdr.loc[keep]
            tmp["comparison"] = comparison
            tmp["region"] = region
            tmp["AS_type"] = ast
            tmp["gene_name"] = tmp.geneSymbol.astype("string").str.strip().str.replace('"', '', regex=False)
            tmp["EventID"] = ast + "|" + tmp[["GeneID"] + coordinates].astype(str).agg("|".join, axis=1)
            event_tables.append(tmp[["region", "comparison", "gene_name", "EventID", "AS_type", "dPSI", "FDR_num"]])

# Save results
as_events = pd.concat(event_tables, ignore_index=True)
as_events["abs_dPSI"] = as_events.dPSI.abs()
as_events.to_csv(OUT_DIR / "AS_testable_events_read10.tsv.gz", sep="\t", index=False)
print(f"Gene records: {len(gene):,}; testable AS events: {len(as_events):,}")
