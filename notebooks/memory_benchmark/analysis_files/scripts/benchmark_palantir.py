#!/usr/bin/env python3

# this requires terminal states from CellRank's GPCCA (see the time benchmarks)


import scanpy as sc
import palantir
import numpy as np
import pandas as pd
import pickle
import os
import traceback

from anndata import AnnData
from typing import Optional, List

from math import ceil


def _clean_orig_names(series):
    return np.array(list(map(np.array, series.str.split(":"))))[:, 1]


def _add_annotations(adata: AnnData, annot: pd.DataFrame) -> None:
    adata.obs["genes_cleaned"] = _clean_orig_names(adata.obs.index)
    assert len(set(adata.obs["genes_cleaned"])) == len(adata)

    annot["genes_cleaned"] = np.array(list(map(np.array, annot.index.str.split("_"))))[
        :, 2
    ]
    annot["genes_cleaned"] = annot["genes_cleaned"].str.replace("-", "x-")

    tmp = adata.obs.merge(annot, how="left", on="genes_cleaned")
    tmp.drop_duplicates("genes_cleaned", inplace=True)
    tmp.set_index("genes_cleaned", drop=True, inplace=True)

    adata.obs = tmp
    adata.obs["Reprogramming Day"] = adata.obs["Reprogramming Day"].astype("category")


def _select_root_cell(adata: AnnData) -> str:
    obs = adata.obs["Reprogramming Day"]
    min_val = np.nanmin(obs.cat.categories)

    return obs[obs == min_val].index[0]


def _load_cellrank_final_states(adata: AnnData, data) -> Optional[list]:
    try:
        # old name: main_states
        index = _clean_orig_names(data["terminal_states"].index)
        valid_ixs = np.isin(index, adata.obs.index)

        x = data["lin_probs"][valid_ixs, :]
        x = pd.DataFrame(x, index=index[valid_ixs])

        if len(index) < 3:
            return None

        ixs = []
        for lin in range(x.shape[1]):
            y = x[~np.isin(x.index, ixs)]

            assert len(y) + len(ixs) == x.shape[0], "Sanity check failed"

            ix = np.argmax(y.values[:, lin])
            ixs.append(y.index[ix])

        return ixs
    except Exception as e:
        print(f"Unexpected error: `{e}`.")
        raise e


def _palantir_preprocess(adata: AnnData):
    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="cell_ranger", n_top_genes=1500)

    print("Running PCA")
    n_comps = 300
    sc.pp.pca(adata, use_highly_variable=True, n_comps=n_comps)

    print("Diff maps")
    dm_res = palantir.utils.run_diffusion_maps(
        pd.DataFrame(adata.obsm["X_pca"][:, :n_comps], index=adata.obs_names)
    )

    print("MS space")
    ms_data = palantir.utils.determine_multiscale_space(dm_res)

    return ms_data


def _benchmark_palantir(
    bdata: AnnData,
    size: int,
    col: int,
    annot: pd.DataFrame,
    fs_data: pd.DataFrame,
    n_jobs: int = 32,
) -> Optional[List[float]]:
    from utils import benchmark

    run_palantir = benchmark(palantir.core.run_palantir)
    res = None

    try:
        print(f"Subsetting data to `{size}`, split `{col}`.")
        _add_annotations(bdata, annot)

        assert bdata.n_obs == size

        root_cell = _select_root_cell(bdata)
        final_states = _load_cellrank_final_states(bdata, fs_data)

        if final_states is None:
            print("No final states found, skipping")
            return None
        if root_cell in final_states:
            print("Root cell is in final states, skipping")
            return None

        print("Preprocessing")
        ms_data = _palantir_preprocess(bdata)

        print(
            f"Running with CellRank terminal states `root_cell={root_cell}` and "
            f"`final_states={final_states}`"
        )
        res, _ = run_palantir(
            ms_data,
            root_cell,
            terminal_states=final_states,
            knn=30,
            num_waypoints=int(ceil(size * 0.15)),
            n_jobs=n_jobs,
            scale_components=False,
            use_early_cell_as_start=True,
        )
    except Exception as e:
        print(
            f"Unable to run `Palantir` with size `{size}` on split `{col}`. Reason: `{e}`."
        )
        print(traceback.format_exc())

    return res


def benchmark_palantir(adata, size: int, col: int, n_jobs: int = 32) -> None:
    from utils import PROFILER_ROOT, PROFILER_ROOT_1_CORE, DATA_DIR

    path = PROFILER_ROOT_1_CORE if n_jobs == 1 else PROFILER_ROOT
    path = path / "palantir" / f"{size}_{col}.pickle"

    if not os.path.isdir(path.parent):
        os.makedirs(path.parent, exist_ok=True)

    annot = pd.read_csv(
        DATA_DIR / "morris_data" / "annotations" / "supp_table_4.csv",
        index_col=0,
        header=2,
    )
    with open(
        DATA_DIR / "benchmarking" / "runtime_analysis" / "gpcca.pickle", "rb"
    ) as fin:
        fs_data = pickle.load(fin)[size][str(col)]  # final states data

    res = _benchmark_palantir(
        adata, size, col, annot=annot, fs_data=fs_data, n_jobs=n_jobs
    )

    if res is not None:
        with open(path, "wb") as fout:
            pickle.dump(res, fout)
