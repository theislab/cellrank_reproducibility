from anndata import AnnData
from pandas.api.types import is_categorical_dtype
from typing import Optional, Union, Iterable, Dict, TypeVar, Tuple, Any
from itertools import combinations
from collections import defaultdict
import matplotlib.cm as cm
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd

from cellrank.tl import Lineage
from cellrank.tl.estimators import BaseEstimator
from scipy.stats import wilcoxon

import cellrank as cr
import scvelo as scv
import scanpy as sc


Cache = TypeVar("Cache")


def _run_cellrank(
    adata: AnnData,
    metast: Optional[int],
    finst: Union[int, Iterable[str]],
    backward: bool = False,
    show_figures: bool = True,
    save_figures: bool = True,
    weight_connectivities: float = 0.2,
    mode: str = "deterministic",
    str_pars: Optional[str] = None,
):
    """Utility function to run CellRank for the forward/backward process."""

    # check whether the types make sense
    if metast is not None and not isinstance(metast, int):
        raise ValueError(f"`metast` must be `int` or `None`, found {type(metast)}")

    if not isinstance(finst, (int, list, pd.Series)) and finst is not None:
        raise ValueError(
            f"`finst` must be `int` or `list` or None, found {type(finst)}"
        )

    if metast is None and not isinstance(finst, (list, pd.Series)):
        raise ValueError("When `metast` is None, `finst` must be a list of str")

    dir_key = "backward" if backward else "forward"

    # initialise kernel object
    kernel = cr.tl.transition_matrix(
        adata,
        backward=backward,
        mode=mode,
        n_jobs=4,
        show_progress_bar=False,
        softmax_scale=None,
        weight_connectivities=weight_connectivities,
    )
    g = cr.tl.estimators.GPCCA(kernel)

    # compute eigendecomposition
    g.compute_schur(which="LM", n_components=20)
    if show_figures:
        figure_params = {"title": f"spectrum {dir_key}", "dpi": 100}
        if save_figures:
            figure_params["save"] = f"spectrum_{dir_key}_{str_pars}.pdf"
        g.plot_spectrum(**figure_params)

    # comptue metastable states. Depending on type(metast), we either compute a fixed number of iterate
    # until the desired states are present
    if isinstance(metast, int):
        n_states = metast
    elif isinstance(finst, pd.Series):
        n_states = len(finst.cat.categories)
    elif metast is None:
        n_states = len(finst)
    else:
        raise RuntimeError("Unable to set `n_states`.")

    while True:
        try:
            g.compute_macrostates(cluster_key="clusters", n_states=n_states)
        except ValueError:
            print(
                f"WARNING: Using one more metastable states than requested, dir={dir_key}"
            )
            n_states += 1
            continue
            # g.compute_metastable_states(cluster_key='clusters', n_states=n_states)

        # check whether all required states are contained in the current set of metastable states
        if metast is None:
            mainst_updated = []
            current_metast = g.macrostates.cat.categories

            for target in (
                finst.cat.categories if is_categorical_dtype(finst) else finst
            ):
                mask = np.array([key.startswith(target) for key in current_metast])
                if not mask.any():
                    print(
                        f"WARNING: Using one more metastable state than requested, dir={dir_key}"
                    )
                    n_states += 1
                    break
                maybe_combined_mainst = ", ".join(current_metast[mask])
                mainst_updated.append(maybe_combined_mainst)
            else:
                break
        else:
            mainst_updated = finst
            break

    # if mainst. is int, we compute the top n states. If it's a list, we set them to these elements
    if is_categorical_dtype(finst):
        g.set_terminal_states(finst)
    elif isinstance(mainst_updated, int):
        g.compute_terminal_states(method="top_n", n_states=mainst_updated)
    else:
        g.set_terminal_states_from_macrostates(mainst_updated)

    g.compute_absorption_probabilities(
        n_jobs=4, backend="loky", show_progress_bar=False, use_petsc=True
    )

    # plot coarse grained transition matrix, if it has been computed
    if g.coarse_T is not None and show_figures:
        g.plot_coarse_T(title=f"{dir_key} T coarse")

        # plot the main states
    if show_figures:
        g.plot_absorption_probabilities(n_cells=30, same_plot=True)
    final_states_names = list(g.terminal_states.cat.categories)

    # in the forwards case, also plot the lineage probs:
    if not backward:

        if show_figures:
            figure_kwargs = {"title": f"cell fates"}
            if save_figures:
                figure_kwargs["save"] = f"lin_{str_pars}"
            g.plot_absorption_probabilities(same_plot=True, **figure_kwargs)

        # Validate lineage probs by comparing with Fev+ subcl.
        if show_figures:
            figure_kwargs = {
                "title": f"average fate probabilities",
                "cluster_key": "clusters_fine",
                "clusters": [
                    "Fev+ Beta",
                    "Fev+ Alpha",
                    "Fev+ Pyy",
                    "Fev+ Delta",
                    "Fev+ Epsilon",
                ],
            }
            if save_figures:
                figure_kwargs["save"] = f"heatmap_{str_pars}.pdf"
            cr.pl.cluster_fates(adata, mode="heatmap", **figure_kwargs)

    return g, final_states_names, n_states


def run_analysis(
    adata: AnnData,
    min_shared_counts: int = 20,
    n_top_genes: int = 2000,
    n_pcs: int = 30,
    n_neighbors: int = 30,
    perc_data: float = 1,
    weight_connectivities: float = 0.2,
    mode: str = "deterministic",
    n_metast_fwd: Optional[int] = 3,
    n_metast_bwd: Optional[int] = 1,
    finst_fwd: Union[int, Iterable[str]] = 3,
    finst_bwd: Union[int, Iterable[str]] = 1,
    seed: int = 0,
    show_figures: bool = False,
    save_figures: bool = True,
    force: bool = False,
    return_objects: bool = True,
    cache_object: Cache = None,
):
    """
    Utility function for robustness analysis.

    Parameters
    ----------
    adata
        Annotated data object.
    min_shared_counts
        Gene filtering threshold for scvelo. Used to identify genes that have enough
        counts in u and s to estimate velocity.
    n_top_genes
        Gene filtering threshold for scvelo. Used to select the number of highly variable genes.
    n_pcs
        Number of principal components to use.
    n_neighbors
        Neighbors for KNN construction.
    perc_data
        If smaller than 1, subsample cell numbers to random sample of given size.
    weight_connectivities
        Weight for :func:`cellrank.tl.transition_matrix`.
    mode
        Velocity kernel mode.
    n_metast_fwd, n_metast_bwd
        If int, compute that number of states. If None, compute enough so that the `mainst` are included.
    finst_fwd, finst_bwd
        If int = x, restrict metastates to most likely x. If list of str, restrict to these. If None, just use the
        metastable states.
    seed
        Seed for random subsampling of the data.
    show_figures
        Show figures like metastable states, spectrum, coarse grained transition matrix etc.
    save_figures
        If True, save the produced figures. Note: If `show_figures = False`, we will also not save them.
    force
        Forced recomputation of dynamical model parameters.
    return_objects
        Wheter to return g_fwd and g_bwd.
    cache_object
        Scachepy cache object.

    Returns
    -------
    results
        dict containing metrics describing the results.
    params
        dict containing the input parameters.

        Also, produces plots of spectra, states and lineages
    """

    # make a copy of the AnnData object
    adata_comp = adata.copy()

    # create a string for saving figures
    str_pars = (
        f"knn_{n_neighbors}_npcs_{n_pcs}_percd_{perc_data}_seed_{seed}_"
        f"nhvg_{n_top_genes}_minsc_{min_shared_counts}_nstat_{n_metast_fwd}_velo_uncert_wc_{weight_connectivities}"
    )
    str_cache = f"knn_{n_neighbors}_npcs_{n_pcs}_percd_{perc_data}_nhvg_{n_top_genes}_minsc_{min_shared_counts}"
    print(f"Running the analysis with `{str_pars}`")

    # potentially downsmaple cell number
    if perc_data < 1:
        print("WARNING: Subsampling the number of cells")
        str_cache += f"_seed_{seed}"
        n_cells = int(perc_data * adata_comp.n_obs)
        np.random.seed(seed)
        bcs = np.random.choice(adata_comp.obs_names, size=n_cells, replace=False)
        if is_categorical_dtype(finst_fwd):
            adata_comp.obs["tmp"] = finst_fwd.values
            adata_comp.obs["tmp"] = adata_comp.obs["tmp"].astype("category")
        print(f"Downsampling to {n_cells} cells")
        adata_comp = adata_comp[bcs].copy()
        finst_fwd = adata_comp.obs["tmp"].astype("category")

    # ------------------------------------------------------
    # pre-processing and scvelo

    # filter, normalize, log transform
    scv.pp.filter_genes(adata_comp, min_shared_counts=min_shared_counts)
    scv.pp.normalize_per_cell(adata_comp)
    scv.pp.filter_genes_dispersion(adata_comp, n_top_genes=n_top_genes)
    scv.pp.log1p(adata_comp)

    # compute pca, knn graph and scvelo's moments
    sc.tl.pca(adata_comp, n_comps=n_pcs)
    sc.pp.neighbors(adata_comp, n_neighbors=n_neighbors, n_pcs=n_pcs)
    scv.pp.moments(adata_comp, n_pcs=n_pcs, n_neighbors=n_neighbors)

    # compute/load from cache the dyn. model params and compute velocities
    if cache_object is not None:
        cache_object.tl.recover_dynamics(
            adata_comp, fname=f"pars_{str_cache}", force=force
        )
    else:
        scv.tl.recover_dynamics(adata_comp)
    scv.tl.velocity(adata_comp, mode="dynamical", min_r2=None)
    scv.tl.velocity_graph(adata_comp)

    # plot streamline plot
    if show_figures:
        figure_params = {
            "color": ["clusters"],
            "title": f"scvelo velocities",
            "legend_loc": "right",
        }
        if save_figures:
            figure_params["save"] = f"velocities_{str_pars}.png"
        scv.pl.velocity_embedding_stream(adata_comp, **figure_params)

    # ------------------------------------------------------
    # cellrank

    # run the forward and backward processes
    g_fwd, terminal_states_names, n_metast_fwd = _run_cellrank(
        adata_comp,
        backward=False,
        mode=mode,
        metast=n_metast_fwd,
        finst=finst_fwd,
        save_figures=save_figures,
        str_pars=str_pars,
        show_figures=show_figures,
        weight_connectivities=weight_connectivities,
    )
    g_bwd, initial_states_names, n_metast_bwd = _run_cellrank(
        adata_comp,
        backward=True,
        mode=mode,
        metast=n_metast_bwd,
        finst=finst_bwd,
        save_figures=save_figures,
        str_pars=str_pars,
        show_figures=show_figures,
        weight_connectivities=weight_connectivities,
    )

    # aggregate root and final states
    # g_fwd.set_final_states(final_states_names, redistribute=False)
    _create_root_final_annotations(adata_comp)
    if show_figures:
        figure_params = {
            "title": f"initial and terminal states",
            "legend_loc": "right",
            "size": 100,
        }
        if save_figures:
            figure_params["save"] = f"initial_terminal_{str_pars}.pdf"
        scv.pl.scatter(adata_comp, color="initial_terminal", **figure_params)

    # ------------------------------------------------------
    # metrics

    e_gap_fwd, e_gap_bwd = (
        g_fwd.eigendecomposition["eigengap"],
        g_bwd.eigendecomposition["eigengap"],
    )
    fev_mask = adata_comp.obs["clusters"] == "Fev+"
    if "rest" in g_fwd.absorption_probabilities.names:
        fev_data_meta = g_fwd.macrostates_memberships[fev_mask, :-1].X
        fev_data = g_fwd.absorption_probabilities[fev_mask, :-1].X
        # fev_fates = np.mean(g_fwd.absorption_probabilities[fev_mask, :-1].X, axis=0)
    else:
        fev_data_meta = g_fwd.macrostates_memberships[fev_mask].X
        fev_data = g_fwd.absorption_probabilities[fev_mask].X
        # fev_fates = np.mean(g_fwd.absorption_probabilities[fev_mask].X, axis=0)
    fev_fates_mean = np.mean(fev_data, axis=0)
    fev_fates_error = np.std(fev_data, axis=0) / np.sqrt(np.sum(fev_mask))

    fev_fates_meta_mean = np.mean(fev_data_meta, axis=0)
    fev_fates_meta_error = np.std(fev_data_meta, axis=0) / np.sqrt(np.sum(fev_mask))

    # construct for dicts and combine them
    dict_1 = {"egap_fwd": e_gap_fwd, "egap_bwd": e_gap_bwd}
    dict_2 = {f"final_{i}_fwd": value for i, value in enumerate(terminal_states_names)}
    dict_3 = {f"final_{i}_fwd_p": value for i, value in enumerate(fev_fates_mean)}
    dict_4 = {f"final_{i}_fwd_err": value for i, value in enumerate(fev_fates_error)}
    dict_5 = {f"meta_{i}_fwd_p": value for i, value in enumerate(fev_fates_meta_mean)}
    dict_6 = {
        f"meta_{i}_fwd_err": value for i, value in enumerate(fev_fates_meta_error)
    }
    dict_7 = {f"final_{i}_bwd": value for i, value in enumerate(initial_states_names)}
    results = {**dict_1, **dict_2, **dict_3, **dict_4, **dict_5, **dict_6, **dict_7}
    params = {
        "n_metas_fwd": n_metast_fwd,
        "n_metas_bwd": n_metast_bwd,
        "min_sh_c": min_shared_counts,
        "n_top_g": n_top_genes,
        "n_pcs": n_pcs,
        "n_neigh": n_neighbors,
        "perc_data": perc_data,
        "seed": seed,
    }

    if return_objects:
        return params, results, g_fwd, g_bwd
    else:
        return params, results


def _pairwise_non_aligned_corr(
    o1: Lineage,
    o2: Lineage,
    *,
    key: str,
    names1: Optional[pd.Index] = None,
    names2: Optional[pd.Index] = None,
):
    """Utility function to handle the case of subsampled data."""

    mask_1 = [state.startswith(key) for state in o1.names]
    mask_2 = [state.startswith(key) for state in o2.names]

    assert np.sum(mask_1) == 1, f"`{key}` not (uniquely) in o1"
    assert np.sum(mask_2) == 1, f"`{key}` not (uniquely) in o2"
    o1.names[mask_1] = key
    o2.names[mask_2] = key

    # construct pandas DataFrames from both
    l1 = pd.DataFrame(o1, index=names1, columns=o1.names)
    l2 = pd.DataFrame(o2, index=names2, columns=o2.names)

    # joint the df's, using 'inner' and index-to-index
    joined_df = l1.join(l2, lsuffix="_l", rsuffix="_r", how="inner")

    # select the columns corresponding to the given keys
    x = joined_df[f"{key}_l"]
    y = joined_df[f"{key}_r"]

    return np.corrcoef(x, y)[0, 1]


def plot_correlation_map(
    data: Dict[Union[str, float], Lineage],
    axes=None,
    show_title: bool = True,
    keys: Tuple[str] = ("Alpha", "Beta", "Epsilon", "Delta"),
    tri: bool = True,
    save: Optional[str] = None,
    barcodes_align: bool = True,
    ylabel: Optional[str] = None,
    show_cbar: bool = True,
    ranges: Optional[Tuple[float, float]] = None,
    return_ranges: bool = False,
    subsampling: Optional[Dict[Union[str, float], Lineage]] = None,
    **kwargs: Any,
) -> Union[Tuple[float, float], Dict[str, pd.DataFrame]]:
    """
    Utility for computing correlation matrices of fates and plotting these.

    Params
    --------
    data
        Dict containing lineage objects as values and perturbed parameter settings as keys.
    keys
        List of strings denoting lineage names.
    tri
        If true, only show half of the correlation matrix (it's symmetric).
    save
        Saving plots using plt.
    barcodes_align
        If True, assumes that all lineage objects have been computed for the same cells in the same order.
        If False, assumes that there is a partial overlap and cells appear in random orders.
    subsampling:
        Dict containing observations names when subsampling approach was used.
    **kwargs
        Keyword arguments for :func:`seaborn.heatmap`.

    Returns
    --------
    Plots a set of heatmaps, one per `key` and returns the correlation matrices.
    """

    if axes is None and not return_ranges:
        fig, axes = plt.subplots(
            1,
            len(keys),
            figsize=(len(keys) * 6, 6),
            sharey="row",
            constrained_layout=True,
        )

    corr_mats = []
    values = {}

    for key in keys:

        if barcodes_align:
            X_list = []
            for obj in data.values():
                # need to be a bit careful here because cyclic states can cause 'or' annotations
                mask = np.array([state.startswith(key) for state in obj.names])
                if np.sum(mask) == 1:
                    X_list.append(obj[:, mask].X)
                elif np.sum(mask) > 1:
                    raise ValueError("Lineage keys are not unique")
                elif np.sum(mask) == 0:
                    X_list.append(np.full_like(obj[:, 0].X, np.nan))

            X = np.concatenate(X_list, axis=1)
            corr_mat = np.corrcoef(X, rowvar=False).T
        else:
            # this case is relevant for subsampled datasets
            corr_mat = np.ones((len(data), len(data)))
            for (j, i), (k1, k2) in zip(
                combinations(range(len(data)), 2), combinations(data.keys(), 2)
            ):
                o1, o2 = data[k1], data[k2]
                if subsampling is None:
                    # assuming names are the same
                    names1, names2 = None, None
                else:
                    names1, names2 = subsampling[k1], subsampling[k2]
                corr_mat[i, j] = _pairwise_non_aligned_corr(o1, o2, key=key, names1=names1, names2=names2)
            corr_mat = corr_mat.T

        corr_mats.append(pd.DataFrame(corr_mat[::-1], columns=data.keys()))

    vmin, vmax = np.min(corr_mats), np.max(corr_mats)
    if return_ranges:
        return vmin, vmax
    elif ranges is not None:
        vmin, vmax = ranges
    cmap = cm.get_cmap("Blues").copy()
    cmap.set_bad(color="whitesmoke")

    for i, (key, ax, corr_mat) in enumerate(zip(keys, axes, corr_mats)):
        if tri:
            mask = np.tril(np.ones_like(corr_mat, dtype=bool), -1)[::-1]
        else:
            mask = None

        sns.heatmap(
            corr_mat,
            vmin=vmin,
            vmax=vmax,
            cmap=cmap,
            ax=ax,
            mask=mask,
            cbar=False,
            square=True,
            **kwargs,
        )
        if show_title:
            ax.set_title(key, y=1.05)
        ax.set_yticks([])
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45)

        quant_mask = ~np.tril(np.ones_like(corr_mat, dtype=bool), 0)[::-1]
        values[key] = corr_mat

        assert np.all(
            ~quant_mask[
                np.arange(quant_mask.shape[0])[::-1], np.arange(quant_mask.shape[0])
            ]
        )

        ax.text(
            0.05,
            0.9,
            f"med = {np.median(corr_mat.values[quant_mask]):.2f}",
            horizontalalignment="left",
            verticalalignment="center",
            transform=ax.transAxes,
        )
        ax.text(
            0.05,
            0.8,
            f"min = {np.min(corr_mat.values[quant_mask]):.2f}",
            horizontalalignment="left",
            verticalalignment="center",
            transform=ax.transAxes,
        )

        if i == 0 and ylabel:
            ax.set_ylabel(ylabel, ha="center", labelpad=20)

    if show_cbar:
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes

        norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
        cax = inset_axes(
            ax,
            width="5%",
            loc="right",
            height="100%",
            bbox_to_anchor=(0.15, 0.0, 1, 1),
            bbox_transform=ax.transAxes,
            borderpad=0,
        )

        cticks = np.linspace(vmin, vmax, 5)

        _ = matplotlib.colorbar.ColorbarBase(
            ax=cax,
            norm=norm,
            cmap=cmap,
            label="correlation",
            ticks=cticks,
            format="%.3f",
        )

    if save is not None:
        fig.savefig(save)

    return values


def _calculate_statistics(
    x: Dict[str, pd.DataFrame], y: Dict[str, pd.DataFrame]
) -> defaultdict:
    """
    Parameters
    ----------
    x
        Stochastic mode samples.
    y
        Deterministic mode samples.
    """
    assert set(x.keys()) == set(y.keys()), "Lineages differ."

    res = defaultdict(dict)
    for alternative in ["greater", "two-sided"]:
        for ln in x.keys():
            _x, _y = x[ln].values, y[ln].values
            assert _x.shape == _y.shape, "Shape mismatch."
            # excludes diagonal, which is all 1s
            mask = ~np.tril(np.ones_like(_x, dtype=bool), 0)[::-1]
            _x = np.ravel(_x[mask])
            _y = np.ravel(_y[mask])
            # Mann-Whitney U is recommended by scipy for > 20 samples
            # For Wilcoxon, "auto" means: #observations < 25 ? exact : approx
            res[alternative][ln] = {
                "res": wilcoxon(_x, _y, mode="auto", alternative=alternative),
                "x": _y,
                "y": _y,
            }

    return res


def _create_root_final_annotations(
    adata: AnnData,
    final_key: str = "terminal_states",
    root_key: str = "initial_states",
    final_pref: Optional[str] = "terminal",
    root_pref: Optional[str] = "initial",
    key_added: Optional[str] = "initial_terminal",
) -> None:
    """
    Create categorical annotations of both root and final states.
    This is a utility function for creating a categorical Series object which combines the information about root
    and final states. The Series is written directly to the AnnData object.  This can for example be used to create a
    scatter plot in scvelo.

    Parameters
    ----------
    adata
        AnnData object to write to (`.obs[key_added]`).
    final_key
        Key from `.obs` where final states have been saved.
    root_key
        Key from `.obs` where root states have been saved.
    final_pref, root_pref
        DirPrefix used in the annotations.
    key_added
        Key added to `adata.obs`.
    Returns
    -------
    Nothing, just writes to AnnData.
    """
    from cellrank.tl._utils import _merge_categorical_series
    from cellrank.tl._colors import _create_categorical_colors

    if f"{final_key}_colors" not in adata.uns:
        adata.uns[f"{final_key}_colors"] = _create_categorical_colors(
            len(adata.obs[final_key].cot.categories)
        )

    if f"{root_key}_colors" not in adata.uns:
        adata.uns[f"{root_key}_colors"] = _create_categorical_colors(30)[::-1][
            len(adata.obs[root_key].cat.categories) :
        ]

    # get both Series objects
    cats_final, colors_final = adata.obs[final_key], adata.uns[f"{final_key}_colors"]
    cats_root, colors_root = adata.obs[root_key], adata.uns[f"{root_key}_colors"]

    # merge
    cats_merged, colors_merged = _merge_categorical_series(
        cats_final, cats_root, list(colors_final), list(colors_root)
    )

    # adjust the names
    final_names = cats_final.cat.categories
    final_labels = [
        f"{final_pref if key in final_names else root_pref}: {key}"
        for key in cats_merged.cat.categories
    ]
    cats_merged.cat.rename_categories(final_labels, inplace=True)

    # write to AnnData
    adata.obs[key_added] = cats_merged
    adata.uns[f"{key_added}_colors"] = colors_merged
