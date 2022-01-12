#!/usr/bin/env python3

import os

import velocyto as vcy

import numpy as np
import pickle
from anndata import AnnData


def _initialize_velocyto(adata: AnnData) -> vcy.VelocytoLoom:
    print("Initializing Velocyto object")

    vlm = object.__new__(vcy.VelocytoLoom)
    vlm.ts = adata.obsm["X_pca"][:, :2]
    vlm.Sx_sz = adata.layers["Ms"].astype(np.float64).T  # needs to be float
    vlm.used_delta_t = 1.0
    vlm.delta_S = adata.layers["velocity"].T
    vlm.S = adata.layers["spliced"].T  # used to get n_neighs
    vlm.delta_S = np.ascontiguousarray(vlm.delta_S)

    return vlm


def _benchmark_velocyto(vlm: vcy.VelocytoLoom, *, n_jobs: int = 32) -> dict:
    from utils import benchmark

    estimate_transition_prob = benchmark(vlm.estimate_transition_prob)
    calculate_embedding_shift = benchmark(vlm.calculate_embedding_shift)
    prepare_markov = benchmark(vlm.prepare_markov)
    run_markov = benchmark(vlm.run_markov)

    print("Calculating transition probabilities")
    etp_mem, _ = estimate_transition_prob(  # only works on 2D embedding
        hidim="Sx_sz",
        embed="ts",
        transform="sqrt",
        psc=1,
        n_neighbors=None,
        knn_random=True,
        n_jobs=n_jobs,
    )
    print("Calculating embedding shift")
    ces_mem, _ = calculate_embedding_shift()

    print("Preparing Markov")
    pm_mem, _ = prepare_markov(sigma_D=1, sigma_W=0.5)
    print("Running Markov")
    rm_mem, _ = run_markov()

    return {
        "estimate_transition_probability": etp_mem,
        "calculate_embedding_shift": ces_mem,
        "prepare_markov": pm_mem,
        "run_markov": rm_mem,
    }


def benchmark_velocyto(adata: AnnData, size: int, col: int, n_jobs: int = 32) -> None:
    from utils import PROFILER_ROOT, PROFILER_ROOT_1_CORE

    path = PROFILER_ROOT_1_CORE if n_jobs == 1 else PROFILER_ROOT
    path = path / "velocyto" / f"{size}_{col}.pickle"

    path.parent.mkdir(parents=True, exist_ok=True)

    res = _benchmark_velocyto(_initialize_velocyto(adata), n_jobs=n_jobs)

    if res is not None:
        with open(path, "wb") as fout:
            pickle.dump(res, fout)
