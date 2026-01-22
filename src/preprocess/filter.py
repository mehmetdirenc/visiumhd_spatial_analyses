from __future__ import annotations
from typing import Optional
import numpy as np
import anndata as ad

def filter_bins(
    adata: ad.AnnData,
    min_counts: int = 20,
    min_genes: int = 10,
    max_pct_mito: Optional[float] = None,
    in_tissue_only: bool = True,
) -> ad.AnnData:
    keep = np.ones(adata.n_obs, dtype=bool)

    if "total_counts" in adata.obs:
        keep &= adata.obs["total_counts"].values >= min_counts
    if "n_genes_by_counts" in adata.obs:
        keep &= adata.obs["n_genes_by_counts"].values >= min_genes

    if max_pct_mito is not None and "pct_mito" in adata.obs:
        keep &= adata.obs["pct_mito"].fillna(0).values <= max_pct_mito

    if in_tissue_only and "in_tissue" in adata.obs:
        keep &= adata.obs["in_tissue"].astype(float).values > 0.5

    return adata[keep].copy()
