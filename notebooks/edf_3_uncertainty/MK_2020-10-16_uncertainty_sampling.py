#!/usr/bin/env python3
# coding: utf-8

# This files produces the sampling part of the suppl. figure involving uncertainty.
# The deterministic and stochastic modes can be found in `ML_2021-10-26_uncertainty.ipynb`.


import cellrank as cr
import scanpy as sc
import scvelo as scv
import numpy as np
import pandas as pd
import pickle
import os
import sys
import gc
import joblib as jl

from time import time
from anndata import AnnData


sys.path.insert(0, "../../")
from paths import DATA_DIR, FIG_DIR

ROOT = DATA_DIR / "benchmarking" / "uncertainty"
n_jobs = 16
n_samples = 50_000

cr.logging.print_versions()

scv.settings.set_figure_params(
    "scvelo", dpi_save=400, dpi=80, transparent=True, fontsize=20, color_map="viridis"
)


if os.path.isfile(ROOT / "adata_preprocessed.h5ad") and os.path.isfile(
    ROOT / "terminal_states.csv"
):
    print(f"Loading data from `{ROOT}`")
    adata = cr.read(ROOT / "adata_preprocessed.h5ad")
    terminal_states = pd.read_csv(ROOT / "terminal_states.csv", index_col=0)["0"]
    terminal_states = terminal_states.astype("category")
else:
    print("Preprocessing")
    adata = cr.datasets.pancreas(DATA_DIR / "pancreas" / "pancreas.h5ad")

    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000, log=True)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    scv.tl.recover_dynamics(adata)

    scv.tl.velocity(adata, mode="dynamical")
    scv.tl.velocity_graph(adata)

    k = cr.tl.transition_matrix(
        adata,
        weight_connectivities=0.2,
        mode="stochastic",
        n_jobs=n_jobs,
        softmax_scale=None,
        show_progress_bar=False,
    )

    g = cr.tl.estimators.GPCCA(k)
    g.compute_schur(20)
    g.compute_macrostates(9, cluster_key="clusters")
    g.set_terminal_states_from_macrostates(["Alpha", "Beta", "Epsilon", "Delta"])

    sc.write(ROOT / "adata_preprocessed.h5ad", adata)
    g.terminal_states.to_csv(ROOT / "terminal_states.csv")
    terminal_states = g.terminal_states


def compute_abs_probs(
    ixs: np.ndarray,
    adata: AnnData,
    terminal_states: pd.Series,
    c: cr.tl.kernels.ConnectivityKernel,
):
    res = []

    for i in ixs:
        try:
            conn = c.copy()
            vk = cr.tl.kernels.VelocityKernel(adata.copy()).compute_transition_matrix(
                seed=i,
                n_jobs=2,
                mode="sampling",
                show_progress_bar=False,
                backend="loky",
                softmax_scale=None,
            )
            k = (0.8 * vk + 0.2 * conn).compute_transition_matrix()

            g = cr.tl.estimators.GPCCA(k)
            g.set_terminal_states(terminal_states)
            g.compute_absorption_probabilities(
                use_petsc=True, n_jobs=2, show_progress_bar=False
            )

            res.append(g.absorption_probabilities.copy())

            del g
            gc.collect()
        except Exception:
            pass

    return res


print("Propagating uncertainty to absorption probabilities")
t = time()
c = cr.tl.kernels.ConnectivityKernel(adata).compute_transition_matrix()
res = jl.Parallel(n_jobs=n_jobs, backend="loky")(
    jl.delayed(compute_abs_probs)(ixs, adata, terminal_states, c)
    for ixs in np.array_split(np.arange(n_samples), n_jobs)
)
print(time() - t)
ln = res[0][0]
res = np.concatenate(res)
print(f"Total: {len(res)} / {n_samples}")

print("Saving")
with open(ROOT / "uncertainty.pickle"), "wb" as fout:
    pickle.dump(res, fout)

mean = cr.tl.Lineage(res.mean(0), names=ln.names, colors=ln.colors)
var = cr.tl.Lineage(res.var(0), names=ln.names, colors=ln.colors)

print("Plotting")
clusters = ["Fev+ Beta", "Fev+ Alpha", "Fev+ Pyy", "Fev+ Delta", "Fev+ Epsilon"]
adata.obsm["to_terminal_states"] = mean

figure_kwargs = {
    "title": f"fate probabilities mean",
    "lineages": ["Alpha", "Beta", "Epsilon", "Delta"],
    "cluster_key": "clusters_fine",
    "clusters": clusters,
    "save": FIG_DIR / "suppl_fig_uncertainty" / "mean.pdf",
    "figsize": (5, 3),
}
cr.pl.cluster_fates(adata, mode="heatmap", **figure_kwargs)


clusters = ["Fev+ Beta", "Fev+ Alpha", "Fev+ Pyy", "Fev+ Delta", "Fev+ Epsilon"]
adata.obsm["to_terminal_states"] = var

figure_kwargs = {
    "title": f"fate probabilities var",
    "lineages": ["Alpha", "Beta", "Epsilon", "Delta"],
    "cluster_key": "clusters_fine",
    "clusters": clusters,
    "save": FIG_DIR / "suppl_fig_uncertainty" / "variance.pdf",
    "figsize": (5, 3),
}
cr.pl.cluster_fates(adata, mode="heatmap", **figure_kwargs)
