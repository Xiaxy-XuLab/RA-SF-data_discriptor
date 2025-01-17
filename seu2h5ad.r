#!/usr/bin/env Rscript
library(SeuratDisk)
library(Seurat)
library(tools)


## load data
args = commandArgs(trailingOnly=TRUE)

inx <- args[1]
sample_index <- file_path_sans_ext(inx)
dat <- readRDS(inx)
tmp_fl <- paste0(sample_index, ".", "h5Seurat")

SaveH5Seurat(dat, tmp_fl)
Convert(paste0(tmp_fl), dest = "h5ad")
system(paste0("rm ", tmp_fl))
