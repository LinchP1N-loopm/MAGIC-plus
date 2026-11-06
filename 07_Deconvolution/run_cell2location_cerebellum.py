#!/usr/bin/env python3
"""Run the complete cell2location workflow for G24-2 cerebellum."""

from __future__ import annotations

import argparse
import gc
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import torch
import cell2location
from cell2location.models import RegressionModel
from cell2location.utils.filtering import filter_genes


# -----------------------------------------------------------------------------
# Fixed project configuration
# -----------------------------------------------------------------------------
BASE = Path("/path/to/project")
REF_PATH = BASE / (
    "ref_datas/GSE165371_cerebellum/"
    "mouse_cerebellum_adult_reference_balanced10k.h5ad"
)
ST_PATH = BASE / "20260210_2/0423_repaired_data/G24-2_repaired.h5ad"
RESULT_DIR = BASE / (
    "20260427_AS_analysis/Fig4_cell2location/results/"
    "cerebellum/G24-2_Cerebellum"
)
REF_DIR = RESULT_DIR.parent / "reference_signatures"
REF_MODEL_DIR = REF_DIR / "model"
ST_MODEL_DIR = RESULT_DIR / "model"
INF_AVER_PATH = REF_DIR / "inf_aver.csv"
OUTPUT_H5AD = RESULT_DIR / "sp_cell2location.h5ad"

LABELS_KEY = "cell_type"
BATCH_KEY = "library"
N_CELLS_PER_LOCATION = 1
GENE_COLUMNS = ["Symbol", "gene_name", "feature_name", "gene_symbol"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref-epochs", type=int, default=250)
    parser.add_argument("--spatial-epochs", type=int, default=30_000)
    parser.add_argument("--spatial-batch-size", type=int, default=1000)
    parser.add_argument("--posterior-samples", type=int, default=100)
    parser.add_argument("--export-batch-size", type=int, default=500)
    parser.add_argument("--force-reference", action="store_true")
    parser.add_argument("--force-spatial", action="store_true")
    return parser.parse_args()


def log(message: object = "") -> None:
    print(message, flush=True)


def model_exists(path: Path) -> bool:
    return (path / "model.pt").exists()


def gene_candidates(adata):
    result = {"var_names": pd.Index(adata.var_names.astype(str))}
    for column in GENE_COLUMNS:
        if column in adata.var.columns:
            values = adata.var[column].astype("string").fillna("").astype(str)
            result[column] = pd.Index(values)
    return result


def apply_gene_names(adata, names):
    names = pd.Series(np.asarray(names), index=np.arange(adata.n_vars)).astype(str)
    valid = ~names.isin(["", "nan", "None", "<NA>"])
    adata = adata[:, valid.to_numpy()].copy()
    adata.var_names = pd.Index(names.loc[valid].to_numpy())
    adata = adata[:, ~adata.var_names.duplicated()].copy()
    adata = adata[:, ~adata.var_names.str.upper().str.startswith("MT-")].copy()
    return adata


def align_reference_and_spatial(adata_sc, adata_st):
    sc_candidates = gene_candidates(adata_sc)
    st_candidates = gene_candidates(adata_st)
    best = None
    for sc_source, sc_names in sc_candidates.items():
        sc_set = set(sc_names) - {"", "nan", "None"}
        for st_source, st_names in st_candidates.items():
            overlap = len(sc_set & (set(st_names) - {"", "nan", "None"}))
            if best is None or overlap > best[0]:
                best = (overlap, sc_source, st_source)
    overlap, sc_source, st_source = best
    log(
        f"Gene namespace: reference={sc_source}, spatial={st_source}, "
        f"overlap={overlap:,}"
    )
    return (
        apply_gene_names(adata_sc, sc_candidates[sc_source]),
        apply_gene_names(adata_st, st_candidates[st_source]),
    )


def align_spatial_to_signature(adata_st, genes):
    candidates = gene_candidates(adata_st)
    gene_set = set(map(str, genes))
    best = max(
        (
            len(gene_set & (set(names) - {"", "nan", "None"})),
            source,
        )
        for source, names in candidates.items()
    )
    overlap, source = best
    log(f"Spatial gene namespace={source}, signature overlap={overlap:,}")
    return apply_gene_names(adata_st, candidates[source])


def prepare_counts(adata, name: str):
    if sp.issparse(adata.X):
        adata.X = adata.X.tocsr()
        values = adata.X.data[:100_000]
    else:
        values = np.asarray(adata.X).ravel()[:100_000]
        adata.X = sp.csr_matrix(adata.X)
    is_integer = np.allclose(values, np.round(values))
    log(f"{name}: shape={adata.shape}, integer_counts={is_integer}")
    if not is_integer:
        raise ValueError(f"{name}.X is not integer raw counts")
    return adata


def extract_signature(adata_ref):
    factor_names = adata_ref.uns["mod"]["factor_names"]
    columns = [f"means_per_cluster_mu_fg_{x}" for x in factor_names]
    if "means_per_cluster_mu_fg" in adata_ref.varm:
        signature = adata_ref.varm["means_per_cluster_mu_fg"][columns].copy()
    else:
        signature = adata_ref.var[columns].copy()
    signature.columns = factor_names
    return signature


def main() -> None:
    args = parse_args()
    REF_DIR.mkdir(parents=True, exist_ok=True)
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    if not REF_PATH.exists():
        raise FileNotFoundError(REF_PATH)
    if not ST_PATH.exists():
        raise FileNotFoundError(ST_PATH)

    if OUTPUT_H5AD.exists() and not args.force_spatial:
        log(f"Final result already exists: {OUTPUT_H5AD}")
        log("Use --force-spatial only if you really want to rerun it.")
        return

    torch.set_float32_matmul_precision("high")
    accelerator = "gpu" if torch.cuda.is_available() else "cpu"
    log(f"cell2location={cell2location.__version__}")
    log(f"device={accelerator}")
    log(f"reference={REF_PATH}")
    log(f"spatial={ST_PATH}")

    adata_st = sc.read_h5ad(ST_PATH)

    # ------------------------------------------------------------------
    # Reference signature: reuse CSV, resume a saved model, or train once.
    # ------------------------------------------------------------------
    if INF_AVER_PATH.exists() and not args.force_reference:
        log(f"Reusing signature: {INF_AVER_PATH}")
        inf_aver = pd.read_csv(INF_AVER_PATH, index_col=0)
        adata_st = align_spatial_to_signature(adata_st, inf_aver.index)
    else:
        adata_sc = sc.read_h5ad(REF_PATH)
        adata_sc = adata_sc[adata_sc.obs[LABELS_KEY].notna()].copy()
        adata_sc, adata_st = align_reference_and_spatial(adata_sc, adata_st)
        adata_sc = prepare_counts(adata_sc, "Reference")

        selected = filter_genes(
            adata_sc,
            cell_count_cutoff=5,
            cell_percentage_cutoff2=0.03,
            nonz_mean_cutoff=1.12,
        )
        adata_ref = adata_sc[:, selected].copy()
        log(f"Filtered reference={adata_ref.shape}")

        RegressionModel.setup_anndata(
            adata=adata_ref,
            batch_key=BATCH_KEY,
            labels_key=LABELS_KEY,
        )

        if model_exists(REF_MODEL_DIR) and not args.force_reference:
            log(f"Resuming saved reference model: {REF_MODEL_DIR}")
            mod_ref = RegressionModel.load(REF_MODEL_DIR, adata=adata_ref)
        else:
            log("Training reference model")
            mod_ref = RegressionModel(adata_ref)
            mod_ref.train(max_epochs=args.ref_epochs, accelerator=accelerator)
            # Save before posterior sampling so an OOM cannot lose training.
            mod_ref.save(REF_MODEL_DIR, overwrite=True)
            log(f"Reference model saved: {REF_MODEL_DIR}")

        log("Exporting reference posterior with low-memory settings")
        adata_ref = mod_ref.export_posterior(
            adata_ref,
            sample_kwargs={
                "num_samples": args.posterior_samples,
                "batch_size": args.export_batch_size,
            },
            export_slot="mod",
            add_to_varm=["means", "stds", "q05", "q95"],
        )
        inf_aver = extract_signature(adata_ref)
        inf_aver.to_csv(INF_AVER_PATH)
        log(f"Signature saved: {INF_AVER_PATH}, shape={inf_aver.shape}")

        del adata_sc, adata_ref, mod_ref
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

    # ------------------------------------------------------------------
    # Spatial model: reuse a saved model if posterior export failed before.
    # ------------------------------------------------------------------
    adata_st = prepare_counts(adata_st, "Spatial")
    common = np.intersect1d(adata_st.var_names, inf_aver.index)
    log(f"Common genes for deconvolution={len(common):,}")
    if len(common) < 3000:
        raise ValueError("Fewer than 3000 common genes")

    adata_st_c2l = adata_st[:, common].copy()
    inf_aver_sub = inf_aver.loc[common].copy()
    cell2location.models.Cell2location.setup_anndata(
        adata=adata_st_c2l,
        batch_key=None,
    )

    if model_exists(ST_MODEL_DIR) and not args.force_spatial:
        log(f"Resuming saved spatial model: {ST_MODEL_DIR}")
        mod_st = cell2location.models.Cell2location.load(
            ST_MODEL_DIR,
            adata=adata_st_c2l,
        )
    else:
        log("Training spatial model")
        mod_st = cell2location.models.Cell2location(
            adata_st_c2l,
            cell_state_df=inf_aver_sub,
            N_cells_per_location=N_CELLS_PER_LOCATION,
            detection_alpha=20,
            detection_mean_per_sample=False,
        )
        batch_size = args.spatial_batch_size or None
        mod_st.train(
            max_epochs=args.spatial_epochs,
            batch_size=batch_size,
            train_size=1,
            accelerator=accelerator,
        )
        # Save before posterior sampling.
        mod_st.save(ST_MODEL_DIR, overwrite=True)
        log(f"Spatial model saved: {ST_MODEL_DIR}")

    log("Exporting spatial posterior with low-memory settings")
    adata_st_c2l = mod_st.export_posterior(
        adata_st_c2l,
        sample_kwargs={
            "num_samples": args.posterior_samples,
            "batch_size": args.export_batch_size,
        },
    )

    factor_names = list(adata_st_c2l.uns["mod"]["factor_names"])
    adata_st_c2l.obs[factor_names] = pd.DataFrame(
        adata_st_c2l.obsm["q05_cell_abundance_w_sf"],
        index=adata_st_c2l.obs_names,
        columns=factor_names,
    )

    temporary = OUTPUT_H5AD.with_name(OUTPUT_H5AD.stem + ".tmp.h5ad")
    adata_st_c2l.write_h5ad(temporary)
    temporary.replace(OUTPUT_H5AD)
    log(f"Completed: {OUTPUT_H5AD}")
    log(f"Result shape={adata_st_c2l.shape}")
    log(f"Cell types ({len(factor_names)}): {factor_names}")


if __name__ == "__main__":
    main()
