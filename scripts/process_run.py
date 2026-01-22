#!/usr/bin/env python3
from pathlib import Path
import argparse
import pandas as pd

from src.utils.log import setup_logging
from src.io.spaceranger import RunSpec, load_run_as_anndata
from src.qc.metrics import add_basic_qc, qc_summary
from src.preprocess.filter import filter_bins

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs_csv", required=True, type=Path)
    ap.add_argument("--run_id", required=True, type=str)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--min_counts", type=int, default=20)
    ap.add_argument("--min_genes", type=int, default=10)
    ap.add_argument("--max_pct_mito", type=float, default=None)
    ap.add_argument("--in_tissue_only", action="store_true")
    args = ap.parse_args()

    setup_logging()
    args.outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.runs_csv)
    row = df.loc[df["run_id"] == args.run_id]
    if row.empty:
        raise SystemExit(f"run_id not found: {args.run_id}")

    r = row.iloc[0]
    run = RunSpec(
        run_id=str(r["run_id"]),
        sample_id=str(r["sample_id"]),
        condition=str(r.get("condition", "")),
        spaceranger_dir=Path(r["spaceranger_dir"]).resolve(),
        bin_um=int(r["bin_um"]),
    )

    adata = load_run_as_anndata(run)
    add_basic_qc(adata)

    summ = qc_summary(adata)
    pd.DataFrame([summ]).to_csv(args.outdir / f"{run.run_id}.qc.csv", index=False)

    adata_f = filter_bins(
        adata,
        min_counts=args.min_counts,
        min_genes=args.min_genes,
        max_pct_mito=args.max_pct_mito,
        in_tissue_only=args.in_tissue_only,
    )

    out_h5ad = args.outdir / f"{run.run_id}.filtered.h5ad"
    adata_f.write(out_h5ad)

    print(f"Wrote: {out_h5ad}")

if __name__ == "__main__":
    main()
