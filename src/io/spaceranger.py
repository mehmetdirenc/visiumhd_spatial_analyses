# src/visiumhd/io/spaceranger.py

from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple, List

import json
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse

# ---------- patterns that match YOUR Space Ranger outputs ----------
# Example:
# outs/binned_outputs/square_024um/filtered_feature_bc_matrix(.h5)
# outs/binned_outputs/square_024um/spatial/tissue_positions.parquet

MATRIX_PATTERNS = [
    # HD binned outputs - YOUR structure
    "outs/binned_outputs/square_{bin3}um/filtered_feature_bc_matrix",
    "outs/binned_outputs/square_{bin3}um/filtered_feature_bc_matrix.h5",
    "outs/binned_outputs/square_{bin3}um/raw_feature_bc_matrix",
    "outs/binned_outputs/square_{bin3}um/raw_feature_bc_matrix.h5",

    # classic 10x (in case you point to non-HD runs)
    "outs/filtered_feature_bc_matrix",
    "outs/filtered_feature_bc_matrix.h5",
]

SPATIAL_PATTERNS = [
    # HD binned outputs - YOUR structure
    "outs/binned_outputs/square_{bin3}um/spatial",

    # fallback
    "outs/spatial",
]

POSITION_FILES = [
    "tissue_positions.parquet",
    "tissue_positions.csv",
    "tissue_positions_list.csv",
]

SCALEFACTORS = "scalefactors_json.json"


@dataclass(frozen=True)
class RunSpec:
    run_id: str
    sample_id: str
    condition: str
    spaceranger_dir: Path
    bin_um: int


def _bin3(bin_um: int) -> str:
    # 2 -> "002", 8 -> "008", 24 -> "024", 100 -> "100"
    return f"{int(bin_um):03d}"


def _resolve_first(base: Path, patterns: List[str], bin_um: int) -> Optional[Path]:
    b3 = _bin3(bin_um)
    for p in patterns:
        cand = base / p.format(bin=bin_um, bin3=b3)
        if cand.exists():
            return cand
    return None


def _read_10x_any(path: Path) -> ad.AnnData:
    if path.is_dir():
        a = sc.read_10x_mtx(str(path), var_names="gene_symbols", make_unique=True)
    elif path.suffix == ".h5":
        a = sc.read_10x_h5(str(path))
    else:
        raise FileNotFoundError(f"Unsupported 10x matrix path: {path}")

    # ensure sparse csr
    a.X = a.X.tocsr() if sparse.issparse(a.X) else sparse.csr_matrix(a.X)
    return a


def _load_spatial(spatial_dir: Path) -> Dict[str, object]:
    out: Dict[str, object] = {}
    if not spatial_dir.exists():
        return out

    sf = spatial_dir / SCALEFACTORS
    if sf.exists():
        out["scalefactors"] = json.loads(sf.read_text())

    pos_path = None
    for fn in POSITION_FILES:
        cand = spatial_dir / fn
        if cand.exists():
            pos_path = cand
            break

    if pos_path is not None:
        if pos_path.suffix == ".parquet":
            pos = pd.read_parquet(pos_path)
        else:
            pos = pd.read_csv(pos_path, header=None)
            if pos.shape[1] >= 6:
                pos.columns = ["barcode","in_tissue","array_row","array_col","pxl_row","pxl_col"] + \
                              [f"extra_{i}" for i in range(pos.shape[1] - 6)]
        out["positions"] = pos

    return out


def load_run_as_anndata(run: RunSpec) -> ad.AnnData:
    base = run.spaceranger_dir

    matrix_path = _resolve_first(base, MATRIX_PATTERNS, run.bin_um)
    if matrix_path is None:
        raise FileNotFoundError(
            f"[{run.run_id}] Could not find matrix for bin {run.bin_um}um under {base}/outs.\n"
            f"Expected something like outs/binned_outputs/square_{_bin3(run.bin_um)}um/filtered_feature_bc_matrix(.h5)"
        )

    a = _read_10x_any(matrix_path)

    # annotate run/sample metadata
    a.obs["run_id"] = run.run_id
    a.obs["sample_id"] = run.sample_id
    a.obs["condition"] = run.condition
    a.uns["bin_um"] = int(run.bin_um)
    a.uns["spaceranger_dir"] = str(base)

    # spatial
    spatial_dir = _resolve_first(base, SPATIAL_PATTERNS, run.bin_um)
    if spatial_dir is not None:
        sp = _load_spatial(spatial_dir)
        if "positions" in sp:
            pos = sp["positions"].copy()

            # handle both common schemas:
            # - column "barcode" (csv style)
            # - index already barcode (some parquet exports)
            if "barcode" in pos.columns:
                pos = pos.set_index("barcode")
            else:
                # if no barcode column, try index as barcode
                pos.index = pos.index.astype(str)

            common = a.obs_names.intersection(pos.index)
            for col in pos.columns:
                a.obs.loc[common, col] = pos.loc[common, col].values

            if {"pxl_row","pxl_col"}.issubset(a.obs.columns):
                a.obsm["spatial"] = a.obs.loc[:, ["pxl_col","pxl_row"]].to_numpy()

        if "scalefactors" in sp:
            a.uns.setdefault("spatial", {})
            a.uns["spatial"][run.sample_id] = {"scalefactors": sp["scalefactors"]}

    return a
