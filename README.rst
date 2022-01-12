CellRank's reproducibility repository
=====================================
We believe that reproducibility is key and have made it as simple as possible to reproduce our results.
Please either open an `issue <https://github.com/theislab/cellrank/issues/new/choose>`_ or contact as
at `info@cellrank.org <mailto:info@cellrank.org>`_ should you experience difficulties reproducing any result.

Manuscript, Code and Data
----------------------
CellRank is published in `Nature Methods`_ and the software package can be found at our main website, `cellrank.org`_. Raw published data is available from the Gene Expression Omnibus (GEO) under accession codes:

- `pancreas`_: `GSE132188 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132188>`_.
- `lung`_: `GSE141259 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141259>`_.
- `reprogramming`_: `GSE99915 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99915>`_.

Processed data, including spliced and unspliced count abundances, is available on
`figshare <https://doi.org/10.6084/m9.figshare.c.5172299>`_.
To ease reproducibility, our data examples can also be accessed trough CellRank's
`dataset interface <https://cellrank.readthedocs.io/en/stable/api.html#module-cellrank.datasets>`_.

Navigating this repository
--------------------------
We've organized this repository along the categories below. For each item, you can click the link under **nbviewer**
to open the notebook in the browser using `nbviewer <https://nbviewer.jupyter.org/>`_.
There is no 1-1 mapping from figures to notebooks - some notebooks produce panels for several figures, and some figures
contain panels from several notebooks.
The tables we provide here make the connection between figures and notebooks explicit. At the top of each notebook,
we indicate the package versions we use. 

Results
-------

.. csv-table:: Main Figures
   :header: "Figure", "nbviewer", "Notebook Path"

    Fig. 1, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/fig_1_concept/ML_2021-09-21_fig_1_concept.ipynb>`__, `path <notebooks/fig_1_concept/ML_2021-09-21_fig_1_concept.ipynb>`__
    Fig. 2, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/fig_2_pancreas_main/ML_2021-09-21_fig_2_and_3_pancreas_main.ipynb>`__, `path <notebooks/fig_2_pancreas_main/ML_2021-09-21_fig_2_and_3_pancreas_main.ipynb>`__
    Fig. 3, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/fig_2_pancreas_main/ML_2021-09-21_fig_2_and_3_pancreas_main.ipynb>`__, `path <notebooks/fig_2_pancreas_main/ML_2021-09-21_fig_2_and_3_pancreas_main.ipynb>`__
    Fig. 4, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/fig_4_mef_reprogramming/ML_2021-09-23_mef_reprogramming.ipynb>`__, `path <notebooks/fig_4_mef_reprogramming/ML_2021-09-23_mef_reprogramming.ipynb>`__
    Fig. 5 *Palantir*, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/fig_5_benchmarking/palantir/ML_2021-10-26_palantir.ipynb>`__, `path <notebooks/fig_5_benchmarking/palantir/ML_2021-10-26_palantir.ipynb>`__
    Fig. 5 *STEMNET*, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/fig_5_benchmarking/stemnet/ML_2020-10-17_plot_fates_and_trends.ipynb>`__, `path <notebooks/fig_5_benchmarking/stemnet/ML_2020-10-17_plot_fates_and_trends.ipynb>`__
    Fig. 5 *Velocyto*, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/fig_5_benchmarking/velocyto/MK_2020-12-01_velocyto.ipynb>`__, `path <notebooks/fig_5_benchmarking/velocyto/MK_2020-12-01_velocyto.ipynb>`__
    Fig. 5 *FateID*, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/fig_5_benchmarking/fateid/ML_2021-10-26_plot_fate_bias.ipynb>`__, `path <notebooks/fig_5_benchmarking/fateid/ML_2021-10-26_plot_fate_bias.ipynb>`__
    Fig. 6, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/fig_6_lung/ML_2021-09-24_fig_6_lung.ipynb>`__, `path <notebooks/fig_6_lung/ML_2021-09-24_fig_6_lung.ipynb>`__

.. csv-table:: Extended Data Figures
   :header: "Figure", "nbviewer", "Notebook Path"

    Extended Data Fig. 1, NA (toy data), NA (toy data)
    Extended Data Fig. 2, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/suppl_fig_GPCCA/ML_2021-10-26_GPCCA.ipynb>`__, `path <notebooks/suppl_fig_GPCCA/ML_2021-10-26_GPCCA.ipynb>`__
    Extended Data Fig. 3, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/edf_3_uncertainty/ML_2021-10-26_uncertainty.ipynb>`__, `path <notebooks/edf_3_uncertainty/ML_2021-10-26_uncertainty.ipynb>`__
    Extended Data Fig. 4, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/fig_2_pancreas_main/ML_2021-09-21_fig_2_and_3_pancreas_main.ipynb>`__, `path <notebooks/fig_2_pancreas_main/ML_2021-09-21_fig_2_and_3_pancreas_main.ipynb>`__
    Extended Data Fig. 5, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/fig_2_pancreas_main/ML_2021-09-21_fig_2_and_3_pancreas_main.ipynb>`__, `path <notebooks/fig_2_pancreas_main/ML_2021-09-21_fig_2_and_3_pancreas_main.ipynb>`__
    Extended Data Fig. 6, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/edf_6_pancreas_ductal/ML_2021-09-22_pancreas_ductal.ipynb>`__, `path <notebooks/edf_6_pancreas_ductal/ML_2021-09-22_pancreas_ductal.ipynb>`__
    Extended Data Fig. 7, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/fig_2_pancreas_main/ML_2021-09-21_fig_2_and_3_pancreas_main.ipynb>`__, `path <notebooks/fig_2_pancreas_main/ML_2021-09-21_fig_2_and_3_pancreas_main.ipynb>`__
    Extended Data Fig. 8, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/fig_2_pancreas_main/ML_2021-09-21_fig_2_and_3_pancreas_main.ipynb>`__, `path <notebooks/fig_2_pancreas_main/ML_2021-09-21_fig_2_and_3_pancreas_main.ipynb>`__
    Extended Data Fig. 9, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/fig_6_lung/ML_2021-09-24_fig_6_lung.ipynb>`__, `path <notebooks/fig_6_lung/ML_2021-09-24_fig_6_lung.ipynb>`__
    Extended Data Fig. 10, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/fig_6_lung/ML_2021-09-24_fig_6_lung.ipynb>`__, `path <notebooks/fig_6_lung/ML_2021-09-24_fig_6_lung.ipynb>`__

.. csv-table:: Supplementary Figures
   :header: "Figure", "nbviewer", "Notebook Path"

    Supplementary Fig. 1, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/fig_2_pancreas_main/ML_2021-09-21_fig_2_and_3_pancreas_main.ipynb>`__, `path <notebooks/fig_2_pancreas_main/ML_2021-09-21_fig_2_and_3_pancreas_main.ipynb>`__
    Supplementary Fig. 2, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/suppl_fig_robustness/MK_2020-10-16_robustness.ipynb>`__, `path <notebooks/suppl_fig_robustness/MK_2020-10-16_robustness.ipynb>`__
    Supplementary Fig. 3, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/fig_2_pancreas_main/ML_2021-09-21_fig_2_and_3_pancreas_main.ipynb>`__, `path <notebooks/fig_2_pancreas_main/ML_2021-09-21_fig_2_and_3_pancreas_main.ipynb>`__
    Supplementary Fig. 4, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/fig_5_benchmarking/palantir/ML_2021-10-26_palantir.ipynb>`__, `path <notebooks/fig_5_benchmarking/palantir/ML_2021-10-26_palantir.ipynb>`__
    Supplementary Fig. 5, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/fig_2_pancreas_main/ML_2021-09-21_fig_2_and_3_pancreas_main.ipynb>`__, `path <notebooks/fig_2_pancreas_main/ML_2021-09-21_fig_2_and_3_pancreas_main.ipynb>`__
    Supplementary Fig. 6, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/suppl_fig_robustness/MK_2020-10-16_robustness.ipynb>`__, `path <notebooks/suppl_fig_robustness/MK_2020-10-16_robustness.ipynb>`__
    Supplementary Fig. 7, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/suppl_fig_robustness/MK_2020-10-16_robustness.ipynb>`__, `path <notebooks/suppl_fig_robustness/MK_2020-10-16_robustness.ipynb>`__
    Supplementary Fig. 8, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/suppl_fig_robustness/MK_2020-10-16_robustness.ipynb>`__, `path <notebooks/suppl_fig_robustness/MK_2020-10-16_robustness.ipynb>`__
    Supplementary Fig. 9, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/suppl_fig_robustness/MK_2020-10-16_robustness.ipynb>`__, `path <notebooks/suppl_fig_robustness/MK_2020-10-16_robustness.ipynb>`__
    Supplementary Fig. 10, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/fig_2_pancreas_main/ML_2021-09-21_fig_2_and_3_pancreas_main.ipynb>`__, `path <notebooks/fig_2_pancreas_main/ML_2021-09-21_fig_2_and_3_pancreas_main.ipynb>`__
    Supplementary Fig. 11, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/fig_2_pancreas_main/ML_2021-09-21_fig_2_and_3_pancreas_main.ipynb>`__, `path <notebooks/fig_2_pancreas_main/ML_2021-09-21_fig_2_and_3_pancreas_main.ipynb>`__
    Supplementary Fig. 12, `link1 <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/fig_2_pancreas_main/ML_2021-09-21_fig_2_and_3_pancreas_main.ipynb>`__ `link2 <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/fig_5_benchmarking/palantir/ML_2021-10-26_palantir.ipynb>`__, `path1 <notebooks/fig_2_pancreas_main/ML_2021-09-21_fig_2_and_3_pancreas_main.ipynb>`__ `path2 <notebooks/fig_5_benchmarking/palantir/ML_2021-10-26_palantir.ipynb>`__
    Supplementary Fig. 13, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/fig_5_benchmarking/fateid/MK_2020-10-17_plot_trends.ipynb>`__, `path <notebooks/fig_5_benchmarking/fateid/MK_2020-10-17_plot_trends.ipynb>`__
    Supplementary Fig. 14, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/fig_5_benchmarking/fateid/MK_2020-10-17_plot_trends.ipynb>`__, `path <notebooks/fig_5_benchmarking/fateid/MK_2020-10-17_plot_trends.ipynb>`__
    Supplementary Fig. 15, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/fig_6_lung/ML_2021-09-24_fig_6_lung.ipynb>`__, `path <notebooks/fig_6_lung/ML_2021-09-24_fig_6_lung.ipynb>`__
    Supplementary Fig. 16, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/fig_6_lung/ML_2021-09-24_fig_6_lung.ipynb>`__, `path <notebooks/fig_6_lung/ML_2021-09-24_fig_6_lung.ipynb>`__
    Supplementary Fig. 17, NA (microscopy results), NA (microscopy results)

.. csv-table:: Supplementary Tables
   :header: "Table", "nbviewer", "Notebook Path"

    Supplementary Tab. 1, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/compute_time_benchmark/MK_2020-10-16_compute_time_benchmark.ipynb>`__, `path <notebooks/compute_time_benchmark/MK_2020-10-16_compute_time_benchmark.ipynb>`__
    Supplementary Tab. 2, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/memory_benchmark/MK_2020-10-16_memory_benchmark.ipynb>`__, `path <notebooks/memory_benchmark/MK_2020-10-16_memory_benchmark.ipynb>`__
    Supplementary Tab. 3, `link <https://nbviewer.org/github/theislab/cellrank_reproducibility/blob/master/notebooks/memory_benchmark/MK_2020-10-16_memory_benchmark_1_core.ipynb>`__, `path <notebooks/memory_benchmark/MK_2020-10-16_memory_benchmark_1_core.ipynb>`__

.. _Nature Methods: https://www.nature.com/articles/s41592-021-01346-6
.. _cellrank.org: https://cellrank.org
.. _pancreas: https://doi.org/10.1242/dev.173849
.. _lung: https://doi.org/10.1038/s41467-020-17358-3
.. _reprogramming: https://doi.org/10.1038/s41586-018-0744-4
