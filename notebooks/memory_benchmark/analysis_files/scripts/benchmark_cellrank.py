#!/usr/bin/env python3
from math import ceil
from anndata import AnnData

import traceback
import os
import pickle

os.environ["NUMBA_DISABLE_CUDA"] = "1"

import scvelo as scv
import cellrank as cr


def create_kernel(adata: AnnData, n_jobs: int = 32, backward: bool = False):
    vk = cr.tl.kernels.VelocityKernel(
        adata, backward=backward
    ).compute_transition_matrix(
        mode="deterministic",
        show_progress_bar=False,
        n_jobs=1 if n_jobs == 1 else 16,
        softmax_scale=4,
        backend="loky",
    )
    ck = cr.tl.kernels.ConnectivityKernel(
        adata, backward=backward
    ).compute_transition_matrix()

    return (0.8 * vk + 0.2 * ck).compute_transition_matrix()


def _benchmark_cellrank(
    bdata, size: int, col: int, Estimator, n_states: int = 3, n_jobs: int = 32
):
    from utils import benchmark

    res = dict()

    print("#jobs:", n_jobs)
    try:
        print(f"Subsetting data to `{size}`, split `{col}`.")
        assert bdata.n_obs == size

        assert bdata.obsp["distances"].shape == (size, size)
        assert bdata.obsp["connectivities"].shape == (size, size)
        assert bdata.uns["velocity_graph"].shape == (size, size)
        assert bdata.uns["velocity_graph_neg"].shape == (size, size)

        print("Recomputing neighbors")
        scv.pp.neighbors(bdata)

        print("Recomputing velocity graph")
        scv.tl.velocity_graph(
            bdata, mode_neighbors="connectivities", n_recurse_neighbors=0
        )

        print("Computing kernel")
        compute_kernel = benchmark(create_kernel)
        kmem, k = compute_kernel(bdata, n_jobs=n_jobs)
        e = Estimator(k)

        compute_ap = benchmark(e.compute_absorption_probabilities)

        if Estimator is cr.tl.estimators.GPCCA:
            compute_schur = benchmark(e.compute_schur)
            compute_macro = benchmark(e.compute_macrostates)
            compute_ld = benchmark(e.compute_lineage_drivers)

            print("Computing Schur decomposition")
            try:
                decmem, _ = compute_schur(n_components=n_states + 1)
            except Exception:
                decmem, _ = compute_schur(n_components=n_states + 2)

            print("Computing macrostates")
            try:
                macromem, _ = compute_macro(n_states=n_states)
            except Exception:
                macromem, _ = compute_macro(n_states=n_states + 1)

            e.set_terminal_states_from_macrostates(n_cells=int(ceil(size // 100)))

            print("Computing absorption probabilities")
            apmem, _ = compute_ap(
                use_petsc=True,
                time_to_absorption=None,
                solver="gmres",
                show_progress_bar=False,
                n_jobs=n_jobs,
                backend="loky",
            )

            drivemem, _ = compute_ld()

            res["dec_mem"] = decmem
            res["macro_mem"] = macromem
            res["kernel_mem"] = kmem
            res["ap_mem"] = apmem
            res["drive_mem"] = drivemem
        elif Estimator is cr.tl.estimators.CFLARE:
            compute_eig = benchmark(e.compute_eigendecomposition)
            compute_terminal = benchmark(e.compute_terminal_states)

            print("Computing eigendecomposition")
            eigmem = compute_eig(k=n_states + 1)

            print("Computing terminal states")
            termmem = compute_terminal(use=n_states)

            print("Computing absorption probabilities")
            apmem = compute_ap(
                use_petsc=True,
                time_to_absorption=None,
                solver="gmres",
                show_progress_bar=False,
                n_jobs=n_jobs,
                backend="loky",
            )
            res["eig_mem"] = eigmem
            res["term_mem"] = termmem
            res["ap_mem"] = apmem
        else:
            raise RuntimeError(f"Invalid Estimator: `{Estimator}`.")

    except Exception as e:
        print(
            f"Unable to run `{Estimator}` with size `{size}` on split `{col}`. Reason: `{e}`."
        )
        print(traceback.format_exc())
        return res

    return res


def benchmark_gpcca(adata: AnnData, size: int, col: int, n_jobs: int = 32):
    from utils import PROFILER_ROOT, PROFILER_ROOT_1_CORE

    path = PROFILER_ROOT_1_CORE if n_jobs == 1 else PROFILER_ROOT
    path = path / "gpcca" / f"{size}_{col}.pickle"

    if not os.path.isdir(path.parent):
        os.makedirs(path.parent, exist_ok=True)

    res = _benchmark_cellrank(adata, size, col, cr.tl.estimators.GPCCA, n_jobs=n_jobs)

    if res is not None:
        with open(path, "wb") as fout:
            pickle.dump(res, fout)
