#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

OUTDIR = Path("../results/filtered")  # adjust
QC = OUTDIR / "qc_summary_all_runs.csv"

df = pd.read_csv(QC)

# ---- QC plots per sample_id ----
metrics = ["n_obs", "median_counts", "median_genes", "median_pct_mito", "frac_in_tissue"]
for sample_id, sub in df.groupby("sample_id"):
    sub = sub.sort_values("bin_um")
    fig, axes = plt.subplots(len(metrics), 1, figsize=(7, 10), sharex=True)
    fig.suptitle(f"{sample_id} — QC vs bin_um")
    for ax, m in zip(axes, metrics):
        if m not in sub.columns:
            ax.axis("off")
            continue
        ax.plot(sub["bin_um"], sub[m], marker="o")
        ax.set_ylabel(m)
    axes[-1].set_xlabel("bin_um")
    fig.tight_layout()
    fig.savefig(OUTDIR / f"qc_vs_bin_{sample_id}.png", dpi=150)
    plt.close(fig)

print("[ok] wrote QC summary plots")

# ---- Spatial quicklooks (counts + genes + in_tissue) for each run ----
# WARNING: if these are huge, limit to a subset first
for run_id in df["run_id"].tolist():
    fp = OUTDIR / f"{run_id}.filtered.h5ad"
    if not fp.exists():
        print("[missing]", fp)
        continue

    adata = sc.read_h5ad(fp)
    # ensure spatial basis exists
    if "spatial" not in adata.obsm:
        print("[no spatial]", run_id)
        continue

    for key in ["total_counts", "n_genes_by_counts", "in_tissue"]:
        if key not in adata.obs:
            continue
        sc.pl.embedding(
            adata,
            basis="spatial",
            color=key,
            show=False,
            title=f"{run_id} — {key}",
        )
        plt.gcf().set_size_inches(6, 5)
        plt.savefig(OUTDIR / f"spatial_{run_id}_{key}.png", dpi=150)
        plt.close()

print("[ok] wrote spatial quicklooks")
