#!/usr/bin/env python3

import os
import sys
import numpy as np
import argparse
import scanpy as sc
import pandas as pd
import pickle
import gc

from typing import Union, Callable
from functools import partial, wraps
from memory_profiler import memory_usage
from pathlib import Path
from anndata import AnnData

from rpy2.robjects import r
import anndata2ri

from benchmark_palantir import benchmark_palantir
from benchmark_cellrank import benchmark_gpcca, benchmark_cflare
from benchmark_velocyto import benchmark_velocyto


sys.path.insert(0, "../../../")
from paths import DATA_DIR


def benchmark(func: Callable, n_jobs: int = 32) -> Callable:

    @wraps(func)
    def wrapper(*args, **kwargs):
        return usage((func, args, kwargs))

    mp = n_jobs != 1
    usage = partial(memory_usage, interval=0.01, include_children=mp, multiprocess=mp, retval=True)

    return wrapper


def save_adata(adata: AnnData, transpose: bool = False):
    anndata2ri.activate()

    if transpose:
        r.saveRDS(adata.X.T, file="adata_t.rds")
    else:
        r.saveRDS(adata.X, file="adata.rds")
    r.saveRDS(adata.obs_names.values, file="obs_names.rds")
    r.saveRDS(adata.var_names.values, file="var_names.rds")

    anndata2ri.deactivate()


def save_stemnet_cluster_pop(size: int, col: int):
    anndata2ri.activate()

    with open(DATA_DIR / "benchmarking" / "runtime_analysis" / "gpcca.pickle", "rb") as fin:
        data = pickle.load(fin)[size][str(col)]

    # old name: main_states
    cluster_annot = data["terminal_states"]
    clusters = cluster_annot.cat.categories

    df = pd.DataFrame(dict(zip(clusters, [cluster_annot.isin([c]) for c in clusters])))
    r.saveRDS(df, file="cluster_pop.rds")

    anndata2ri.deactivate()


PROFILER_ROOT = DATA_DIR / "benchmarking" / "memory_analysis"
PROFILER_ROOT_1_CORE = DATA_DIR / "benchmarking" / "memory_analysis_1_core"
SPLIT_PATH = DATA_DIR / "morris_data" / "splits"

PREPROCESSED_PATH = DATA_DIR / "morris_data" / "adata_preprocessed.h5ad"
RAW_PATH = DATA_DIR / "morris_data" / "adata.h5ad"


def _get_split(dir: Union[Path, str]) -> dict:
    dfs = {}
    split_root = Path(dir)

    for split in os.listdir(split_root):
        if split.endswith('.csv'):
            size = int(split[:-4].split('_')[1])
            dfs[size] = pd.read_csv(split_root / split, index_col=0)

    return {k: dfs[k] for k in sorted(dfs.keys()) if k != 104679}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("method", type=str, choices=["gpcca", "cflare", "palantir", "stemnet", "fateid", "velocyto"],
                        metavar="METHOD")
    parser.add_argument("size", type=int, choices=[i * 10_000 for i in range(1, 11)],
                        metavar="SIZE")
    parser.add_argument("split", type=int, choices=list(range(10)),
                        metavar="SPLIT")
    parser.add_argument("--n-jobs", type=int, default=32)

    args = parser.parse_args()

    print(f"Method: `{args.method}`, size: `{args.size}`, split: `{args.split}`, jobs: `{args.n_jobs}`.")

    if args.method in ("fateid", "palantir"):
        adata = sc.read(RAW_PATH)
    else:
        adata = sc.read(PREPROCESSED_PATH)

    ixs = _get_split(SPLIT_PATH)[args.size][str(args.split)].values
    assert isinstance(ixs, np.ndarray)
    assert ixs.shape == (args.size, )

    adata = adata[ixs].copy()
    gc.collect()

    if args.method == "gpcca":
        benchmark_gpcca(adata, args.size, args.split, args.n_jobs)
    elif args.method == "cflare":
        benchmark_cflare(adata, args.size, args.split, args.n_jobs)
    elif args.method == "palantir":
        benchmark_palantir(adata, args.size, args.split, args.n_jobs)
    elif args.method == "velocyto":
        benchmark_velocyto(adata, args.size, args.split, args.n_jobs)
    # the below 2 just prepare the data
    elif args.method == "stemnet":
        os.makedirs(PROFILER_ROOT / "stemnet", exist_ok=True)
        save_adata(adata)
        save_stemnet_cluster_pop(args.size, args.split)
    elif args.method == "fateid":
        os.makedirs(PROFILER_ROOT / "fateid", exist_ok=True)
        save_adata(adata, transpose=True)
    else:
        raise NotImplementedError(f"Method `{args.method!r}` is not implemented.")

