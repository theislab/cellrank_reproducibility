#!/usr/bin/env Rscript

library(SparseM)
library(peakRAM)
library(FateID)
library(RaceID)

print(packageVersion("FateID"))
print(packageVersion("RaceID"))

args = commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2)

size <- as.numeric(args[1])
split <- as.numeric(args[2])
profiler_root <- "../../../data/benchmarking/memory_analysis/fateid/"
fname = paste0(profiler_root, paste("fateid", size, split, sep = "_"), ".txt")
drivers <- paste0("../../../data/morris_data/", "lin_drivers/")

select_driver_genes <- function(valid_genes, ld) {
    ld <- ld[ld[[1]] %in% valid_genes,]
    selected_genes <- list()
    
    for (i in c(2, 3, 4)) {
        c <- colnames(ld)[i]
        genes <- unlist(selected_genes)
        y <- ld[!(ld[[1]] %in% genes), ]
        y <- y[order(y[, c], decreasing=TRUE),]
        selected_genes[[c]] <- as.vector(y$X[seq(1, 3)])
    }
    
    selected_genes
}

print("Reading data")
X <- as(readRDS("adata_t.rds"), "dgCMatrix")
rownames(X) <- readRDS("var_names.rds")
colnames(X) <- readRDS("obs_names.rds")

stopifnot(dim(X)[2] == size)
mintotal <- size %/% 1000


print("Preprocessing")

sc <- SCseq(X)
sc <- filterdata(sc, mintotal=mintotal)
x <- getfdata(sc)[sc@cluster$features,]

print("Getting marker genes")

valid_genes <- rownames(x)
ld <- read.csv(paste0(drivers, size, "_", split, ".csv"))
markers <- select_driver_genes(valid_genes, ld)

minnr = minnrh = size %/% 100
n <- size %/% 500

pa <- getPart(x, markers, n=n)
clustering <- pa$part
endpoints <- pa$tar
z <- sc@distances

print("Profiling fateBias")
x <- as.matrix(x)

mem <- peakRAM({
    fb <- fateBias(x, clustering,
                   endpoints, z=z,
                   minnr=minnr, minnrh=minnrh,
                   seed=123, use.dist=FALSE)
})
sink(fname)
print(mem)
sink()
