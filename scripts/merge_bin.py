#!/usr/bin/env python3
from pathlib import Path
import argparse
import pandas as pd
import scanpy as sc

from src.integrate.merge import merge_runs

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs_csv", required=True, type=Path)
    ap.add_argument("--filtered_dir", required=True, type=Path)
    ap.add_argument("--bin_um", required=True, type=int)
    ap.add_argument("--out_h5ad", required=True, type=Path)
    args = ap.parse_args()

    df = pd.read_csv(args.runs_csv)
    df = df[df["bin_um"] == args.bin_um].copy()
    if df.empty:
        raise SystemExit(f"No runs found with bin_um={args.bin_um}")

    adatas = []
    for run_id in df["run_id"].tolist():
        fp = args.filtered_dir / f"{run_id}.filtered.h5ad"
        adatas.append(sc.read_h5ad(fp))

    merged = merge_runs(adatas, join="outer")
    merged.uns["bin_um"] = int(args.bin_um)
    merged.write(args.out_h5ad)
    print(f"Wrote merged: {args.out_h5ad}")

if __name__ == "__main__":
    main()
