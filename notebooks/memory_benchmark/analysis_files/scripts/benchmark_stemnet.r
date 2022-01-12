#!/usr/bin/env Rscript

library(peakRAM)
library(STEMNET)
print(packageVersion("STEMNET"))

args = commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2)

size <- args[1]
split <- args[2]
profiler_root <- "../../../data/benchmarking/memory_analysis/stemnet/"
fname = paste0(profiler_root, paste("stemnet", size, split, sep = "_"), ".txt")

expression <- as.matrix(readRDS("adata.rds"))
rownames(expression) <- readRDS("obs_names.rds")
colnames(expression) <- readRDS("var_names.rds")
cluster_pop <- readRDS("cluster_pop.rds")

pop <- booleanTable2Character(cluster_pop, other_value=NA)

print(paste("Profiling STEMNET:", size, split, sep=" "))

mem <- peakRAM({result <- runSTEMNET(expression, pop)})
sink(fname)
print(mem)
sink()
