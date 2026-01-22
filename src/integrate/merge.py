from __future__ import annotations
from typing import List
import anndata as ad

def merge_runs(adatas: List[ad.AnnData], join: str = "outer") -> ad.AnnData:
    # label already in obs; concat keeps run_id/sample_id columns
    merged = ad.concat(adatas, join=join, axis=0, merge="same")
    return merged
