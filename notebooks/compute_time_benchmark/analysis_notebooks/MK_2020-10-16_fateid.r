#!/usr/bin/env Rscript
# Benchmark FateID runtime
# ---
#
# In this notebook, we benchmark the runtime of FateID's `fateBias` function.

# # Preliminaries

# ## Dependency notebooks

# 1. [../preprocessing_notebooks/MK_2020-10-16_generate_splits.ipynb](../preprocessing_notebooks/MK_2020-10-16_generate_splits.ipynb) - extracts the raw matrix, subsets and splits and
# the genes and cell names as .csv files.
# 2. [utils/MK_2020-10-16_save_drivers.ipynb](utils/MK_2020-10-16_save_drivers.ipynb) - extracts the driver genes for each lineage using CellRanks' GPCCA.

# ## Import packages

library(Matrix)
library(SparseM)
library(FateID)
library(RaceID)
library(R.utils)

print(packageVersion("FateID"))
print(packageVersion("RaceID"))


options("scipen"=10)
sizes <- seq(10000, 100000, 10000)

select_driver_genes <- function(valid_genes, ld) {
    selected_genes = {}
    
    ld <- ld[ld[[1]] %in% valid_genes,]
    selected_genes <- list()
    
    for (i in c(2, 3, 4)) {
        c <- colnames(ld)[i]
        genes <- unlist(selected_genes)
        y <- ld[!(ld[[1]] %in% genes), ]
        y <- y[order(y[, c], decreasing=TRUE),]  # select the driver genes with highest correlation
        selected_genes[[c]] <- as.vector(y$X[seq(1, 3)])
    }
    
    selected_genes
}

root_dir <- "../../../data/morris_data/"

print("Reading data")
mat <- Matrix::readMM(paste(root_dir, "raw.mtx", sep=""))
mat <- as(mat, "dgCMatrix")
dim(mat)

res <- matrix(nrow=10, ncol=length(sizes))
colnames(res) <- sizes
rownames(res) <- seq(1, 10)

split_dir <- paste(root_dir, "splits/", sep="")
drivers <- paste(root_dir, "lin_drivers/", sep="")

genes <- read.csv(paste(root_dir, "annotations/genes.csv", sep=""), col.names=c('gene'), skip=0)[['gene']]
obs <- read.csv(paste(root_dir, "annotations/obs_names.csv", sep=""), col.names=c('obs'), skip=0)[['obs']]

rownames(mat) <- genes
colnames(mat) <- obs


split_rows = sapply(seq(0, 9), function(i) {paste("X", i, sep="")})

for (i in seq_along(sizes)) {
    size <- sizes[i]
    fname <- paste(split_dir, "size_", size, ".csv", sep='')
    print(fname)
    splits <- read.csv(fname, header=TRUE) + 1
    splits <- splits[split_rows]
    
    tryCatch({
        for (j in seq(1, 10)) {
            print("Split", j)
            ixs <- splits[[j]]
            stopifnot(length(ixs) == size)
            X <- mat[,ixs]
            
            stopifnot(dim(X)[2] == size)
            mintotal = size %/% 1000
            
            print("Preprocessing")
            sc <- SCseq(X)
            sc <- filterdata(sc, mintotal=mintotal)
            x <- getfdata(sc)[sc@cluster$features,]
            
            print("Getting marker genes")
            valid_genes <- rownames(x)
            ld <- read.csv(paste(drivers, size, "_", j - 1, ".csv", sep=""))
            markers <- select_driver_genes(valid_genes, ld)
                           
            minnr = minnrh = size %/% 100
            n = size %/% 500
                           
            pa <- getPart(x, markers, n=n)
            clustering <- pa$part
            endpoints <- pa$tar
            z <- sc@distances

            print("Running fate bias")
            runtime  <- withTimeout({{
                x <- as.matrix(x)
                start_time <- Sys.time()
                fb <- fateBias(x, clustering,
                               endpoints, z=z,
                               minnr=minnr, minnrh=minnrh,
                               seed=123, use.dist=FALSE)
                end_time <- difftime(Sys.time(), start_time, units="secs")
                end_time
                }},
                timeout=60 * 60 * 4,  # 4 hours
                onTimeout="silent"
            )
            if (is.null(runtime)) {  # timeout
                 runtime <- NA
            }
            
            res[j, i] = runtime
            write.csv(res, "../../../data/benchmarking/runtime_analysis/fateid.csv")
        }
    },
    error=function(e) {print(e)}
    )
}
