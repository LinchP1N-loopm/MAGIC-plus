#!/usr/bin/env python3
"""Use one Zeng Aging reference to deconvolve six G21 spatial samples."""

from __future__ import annotations

import argparse
import gc
from pathlib import Path

import cell2location
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import torch
from cell2location.models import Cell2location, RegressionModel
from cell2location.utils.filtering import filter_genes


SAMPLES = ["G21-4", "G21-7", "G21-1", "G21-8", "G21-2", "G21-5"]
BASE = Path("/path/to/project")
REF_PATH = BASE / (
    "ref_datas/2025_Nature_Zeng_Aging/"
    "Zeng-Aging-Mouse-cell2location-reference-balanced.h5ad"
)
SPATIAL_DIR = BASE / "20260210_2/pp_data"

# 与旧 reference 的结果彻底分开。
OUTPUT_DIR = BASE / "20260705_Fig3/Fig3_cell2location/Zeng_Aging_reference"
REF_DIR = OUTPUT_DIR / "reference_signatures"
REF_MODEL_DIR = REF_DIR / "model"
SIGNATURE_PATH = REF_DIR / "inf_aver.csv.gz"

CELLTYPE_KEY = "cell_type"
BATCH_KEY = "library_label"
MIN_SHARED_GENES = 3_000
N_CELLS_PER_LOCATION = 1
DETECTION_ALPHA = 20
RANDOM_SEED = 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref-epochs", type=int, default=250)
    parser.add_argument("--spatial-epochs", type=int, default=30_000)
    parser.add_argument("--posterior-samples", type=int, default=100)
    parser.add_argument("--export-batch-size", type=int, default=5_000)
    parser.add_argument("--samples", nargs="+", choices=SAMPLES, default=SAMPLES)
    parser.add_argument("--force-reference", action="store_true")
    parser.add_argument("--force-spatial", action="store_true")
    parser.add_argument("--check-only", action="store_true")
    return parser.parse_args()


def log(message: object = "") -> None:
    print(message, flush=True)


def model_exists(path: Path) -> bool:
    return (path / "model.pt").exists()


def validate_counts(adata, name: str, max_values: int = 500_000) -> None:
    if sp.issparse(adata.X):
        adata.X = adata.X.tocsr()
        values = adata.X.data
    else:
        adata.X = sp.csr_matrix(adata.X)
        values = adata.X.data
    if values.size > max_values:
        rng = np.random.default_rng(RANDOM_SEED)
        values = values[rng.choice(values.size, max_values, replace=False)]
    valid = values.size > 0 and values.min() >= 0 and np.allclose(values, np.rint(values))
    log(f"{name}: shape={adata.shape}, integer_counts={valid}")
    if not valid:
        raise ValueError(f"{name}.X is not non-negative integer raw counts")


def standardise_gene_ids(adata):
    adata.var_names = pd.Index(adata.var_names.astype(str)).str.replace(r"\.\d+$", "", regex=True)
    if not adata.var_names.is_unique:
        duplicated = adata.var_names[adata.var_names.duplicated()].unique()[:10].tolist()
        raise ValueError(f"Duplicated gene IDs after version removal: {duplicated}")
    symbol_col = next((x for x in ["Symbol", "gene_symbol", "gene_name"] if x in adata.var), None)
    if symbol_col is not None:
        keep = ~adata.var[symbol_col].astype(str).str.upper().str.startswith("MT-").to_numpy()
        adata = adata[:, keep].copy()
    return adata


def extract_signature(adata_ref) -> pd.DataFrame:
    factors = list(adata_ref.uns["mod"]["factor_names"])
    columns = [f"means_per_cluster_mu_fg_{x}" for x in factors]
    if "means_per_cluster_mu_fg" in adata_ref.varm:
        signature = adata_ref.varm["means_per_cluster_mu_fg"][columns].copy()
    else:
        signature = adata_ref.var[columns].copy()
    signature.columns = factors
    return signature


def train_or_load_reference(args, accelerator: str) -> pd.DataFrame:
    REF_DIR.mkdir(parents=True, exist_ok=True)
    if SIGNATURE_PATH.exists() and not args.force_reference:
        log(f"Reusing reference signature: {SIGNATURE_PATH}")
        return pd.read_csv(SIGNATURE_PATH, index_col=0)

    log(f"Reading Zeng Aging reference: {REF_PATH}")
    adata_sc = standardise_gene_ids(sc.read_h5ad(REF_PATH))
    validate_counts(adata_sc, "Zeng Aging reference")
    for key in [CELLTYPE_KEY, BATCH_KEY]:
        if key not in adata_sc.obs:
            raise KeyError(f"Reference obs lacks {key!r}")
    if adata_sc.obs[CELLTYPE_KEY].isna().any():
        raise ValueError("Reference contains missing cell_type values")
    log(
        f"Reference cells={adata_sc.n_obs:,}, genes={adata_sc.n_vars:,}, "
        f"cell_types={adata_sc.obs[CELLTYPE_KEY].nunique()}, "
        f"adult={(adata_sc.obs['donor_age_category'].astype(str) == 'adult').sum():,}, "
        f"aged={(adata_sc.obs['donor_age_category'].astype(str) == 'aged').sum():,}"
    )

    selected = filter_genes(
        adata_sc,
        cell_count_cutoff=5,
        cell_percentage_cutoff2=0.03,
        nonz_mean_cutoff=1.12,
    )
    adata_ref = adata_sc[:, selected].copy()
    log(f"Filtered reference={adata_ref.shape}")
    RegressionModel.setup_anndata(adata_ref, batch_key=BATCH_KEY, labels_key=CELLTYPE_KEY)

    if model_exists(REF_MODEL_DIR) and not args.force_reference:
        log(f"Resuming saved reference model: {REF_MODEL_DIR}")
        model = RegressionModel.load(REF_MODEL_DIR, adata=adata_ref)
    else:
        log(f"Training reference model for {args.ref_epochs} epochs")
        model = RegressionModel(adata_ref)
        model.train(max_epochs=args.ref_epochs, accelerator=accelerator)
        model.save(REF_MODEL_DIR, overwrite=True)
        log(f"Reference model saved: {REF_MODEL_DIR}")

    log("Exporting reference posterior")
    adata_ref = model.export_posterior(
        adata_ref,
        sample_kwargs={
            "num_samples": args.posterior_samples,
            "batch_size": min(args.export_batch_size, adata_ref.n_obs),
        },
        export_slot="mod",
        add_to_varm=["means", "stds", "q05", "q95"],
    )
    signature = extract_signature(adata_ref)
    signature.to_csv(SIGNATURE_PATH, compression="gzip")
    log(f"Signature saved: {SIGNATURE_PATH}, shape={signature.shape}")
    del adata_sc, adata_ref, model
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    return signature


def deconvolve_sample(sample: str, signature: pd.DataFrame, args, accelerator: str) -> None:
    input_path = SPATIAL_DIR / f"{sample}_repaired.h5ad"
    sample_dir = OUTPUT_DIR / sample
    model_dir = sample_dir / "model"
    output_path = sample_dir / f"{sample}_Zeng_Aging_cell2location.h5ad"
    sample_dir.mkdir(parents=True, exist_ok=True)

    if output_path.exists() and not args.force_spatial:
        log(f"[{sample}] final output exists, skipping: {output_path}")
        return
    if not input_path.exists():
        raise FileNotFoundError(input_path)

    log(f"[{sample}] reading {input_path}")
    adata = sc.read_h5ad(input_path)
    if "in_tissue" not in adata.obs:
        raise KeyError(f"[{sample}] obs lacks 'in_tissue'")
    adata = adata[adata.obs["in_tissue"].astype(str) == "1"].copy()
    if adata.n_obs == 0:
        raise ValueError(f"[{sample}] no spots remain after in_tissue == 1")
    adata = standardise_gene_ids(adata)
    validate_counts(adata, sample)

    shared = adata.var_names[adata.var_names.isin(signature.index)]
    log(f"[{sample}] in-tissue spots={adata.n_obs:,}; shared genes={len(shared):,}")
    if len(shared) < MIN_SHARED_GENES:
        raise ValueError(f"[{sample}] only {len(shared)} shared genes")
    adata_map = adata[:, shared].copy()
    signature_sub = signature.loc[shared].copy()
    Cell2location.setup_anndata(adata_map, batch_key=None)

    if model_exists(model_dir) and not args.force_spatial:
        log(f"[{sample}] resuming saved spatial model: {model_dir}")
        model = Cell2location.load(model_dir, adata=adata_map)
    else:
        log(f"[{sample}] training spatial model full-batch for {args.spatial_epochs} epochs")
        model = Cell2location(
            adata_map,
            cell_state_df=signature_sub,
            N_cells_per_location=N_CELLS_PER_LOCATION,
            detection_alpha=DETECTION_ALPHA,
            detection_mean_per_sample=False,
        )
        # Full-batch: the complete set of spots for this sample is used per epoch.
        model.train(
            max_epochs=args.spatial_epochs,
            batch_size=None,
            train_size=1,
            accelerator=accelerator,
        )
        model.save(model_dir, overwrite=True)
        log(f"[{sample}] spatial model saved: {model_dir}")

    log(f"[{sample}] exporting posterior")
    adata_map = model.export_posterior(
        adata_map,
        sample_kwargs={
            "num_samples": args.posterior_samples,
            "batch_size": min(args.export_batch_size, adata_map.n_obs),
        },
    )
    factors = list(adata_map.uns["mod"]["factor_names"])
    abundance = adata_map.obsm["q05_cell_abundance_w_sf"]
    if not isinstance(abundance, pd.DataFrame):
        abundance = pd.DataFrame(abundance, index=adata_map.obs_names, columns=factors)
    else:
        if abundance.shape[1] != len(factors):
            raise ValueError("Posterior abundance column count does not match factor_names")
        abundance = abundance.copy()
        abundance.columns = factors
    relative = abundance.div(abundance.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    adata_map.obsm["q05_cell_abundance_w_sf"] = abundance
    adata_map.obsm["q05_relative_abundance"] = relative
    adata_map.obs[factors] = abundance
    adata_map.obs["Cell_type_top"] = abundance.idxmax(axis=1).astype("category")
    adata_map.obs["Cell_type_top_fraction"] = relative.max(axis=1)
    adata_map.uns["deconvolution_reference"] = {
        "name": "Zeng Aging Mouse balanced reference",
        "cell_type_level": "WMB subclass",
        "reference_path": str(REF_PATH),
        "sample": sample,
    }

    summary = pd.DataFrame({
        "q05_abundance_sum": abundance.sum(axis=0),
        "q05_abundance_mean": abundance.mean(axis=0),
        "mean_relative_abundance": relative.mean(axis=0),
    }).sort_values("q05_abundance_sum", ascending=False)
    abundance.to_csv(sample_dir / f"{sample}_cell_abundance_q05.tsv.gz", sep="\t", compression="gzip")
    relative.to_csv(
        sample_dir / f"{sample}_cell_abundance_q05_relative.tsv.gz",
        sep="\t",
        compression="gzip",
    )
    summary.to_csv(sample_dir / f"{sample}_celltype_abundance_summary.tsv", sep="\t")
    temporary = output_path.with_name(output_path.stem + ".tmp.h5ad")
    adata_map.write_h5ad(temporary)
    temporary.replace(output_path)
    log(f"[{sample}] completed: {output_path}")

    del adata, adata_map, model, abundance, relative
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()


def check_inputs(samples: list[str]) -> None:
    if not REF_PATH.exists():
        raise FileNotFoundError(REF_PATH)
    adata_ref = sc.read_h5ad(REF_PATH, backed="r")
    log(
        f"Reference check: shape={adata_ref.shape}, "
        f"cell_types={adata_ref.obs[CELLTYPE_KEY].nunique()}, "
        f"adult={(adata_ref.obs['donor_age_category'].astype(str) == 'adult').sum():,}, "
        f"aged={(adata_ref.obs['donor_age_category'].astype(str) == 'aged').sum():,}"
    )
    reference_genes = adata_ref.var_names.copy()
    adata_ref.file.close()
    for sample in samples:
        path = SPATIAL_DIR / f"{sample}_repaired.h5ad"
        if not path.exists():
            raise FileNotFoundError(path)
        adata = sc.read_h5ad(path, backed="r")
        in_tissue = int((adata.obs["in_tissue"].astype(str) == "1").sum())
        shared = len(adata.var_names.intersection(reference_genes))
        log(f"[{sample}] shape={adata.shape}, in_tissue={in_tissue:,}, shared genes={shared:,}")
        adata.file.close()


def main() -> None:
    args = parse_args()
    if args.force_reference:
        args.force_spatial = True
    np.random.seed(RANDOM_SEED)
    torch.manual_seed(RANDOM_SEED)
    torch.set_float32_matmul_precision("high")
    accelerator = "gpu" if torch.cuda.is_available() else "cpu"
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    log(f"cell2location={cell2location.__version__}")
    log(f"device={accelerator}")
    log(f"reference={REF_PATH}")
    log(f"output={OUTPUT_DIR}")
    if args.check_only:
        check_inputs(args.samples)
        log("CHECK OK")
        return
    signature = train_or_load_reference(args, accelerator)
    for sample in args.samples:
        deconvolve_sample(sample, signature, args, accelerator)
    log("All requested samples completed.")


if __name__ == "__main__":
    main()

