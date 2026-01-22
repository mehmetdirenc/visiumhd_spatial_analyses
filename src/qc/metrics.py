from __future__ import annotations
from typing import Tuple, Dict, Any
import numpy as np
import scanpy as sc
import anndata as ad

def add_basic_qc(adata: ad.AnnData, mito_prefix: Tuple[str, ...] = ("MT-","mt-")) -> None:
    sc.pp.calculate_qc_metrics(adata, inplace=True)

    mito = np.zeros(adata.n_vars, dtype=bool)
    for p in mito_prefix:
        mito |= adata.var_names.str.startswith(p)
    if mito.any():
        mito_counts = adata.X[:, mito].sum(axis=1)
        total = adata.X.sum(axis=1)
        adata.obs["pct_mito"] = 100.0 * (np.asarray(mito_counts).ravel() / (np.asarray(total).ravel() + 1e-9))
    else:
        adata.obs["pct_mito"] = np.nan

def qc_summary(adata: ad.AnnData) -> Dict[str, Any]:
    out = {
        "run_id": adata.obs["run_id"].iloc[0],
        "sample_id": adata.obs["sample_id"].iloc[0],
        "condition": adata.obs["condition"].iloc[0],
        "bin_um": adata.uns.get("bin_um", np.nan),
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "median_counts": float(np.median(adata.obs["total_counts"])) if "total_counts" in adata.obs else np.nan,
        "median_genes": float(np.median(adata.obs["n_genes_by_counts"])) if "n_genes_by_counts" in adata.obs else np.nan,
        "median_pct_mito": float(np.nanmedian(adata.obs["pct_mito"])) if "pct_mito" in adata.obs else np.nan,
    }
    if "in_tissue" in adata.obs:
        out["frac_in_tissue"] = float(np.mean(adata.obs["in_tissue"].astype(float) > 0.5))
    return out
