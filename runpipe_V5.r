#!/usr/bin/env R
# coding utf-8

###################### Package Loading ###########################
pacman::p_load(
  "Seurat", "SeuratWrappers", "ggplot2",
  "batchelor", "harmony",
  "dplyr", "optparse", "reshape2",
  "data.table", "magrittr", "log4r"
)
#############################################################


###################### META ARGVs ###########################
#############################################################
mc.cores <- 10
############################################################

########################## FUN ##############################
#############################################################
.read10x <- function(
    d, sample_index = NULL, num_min_cell = 3,
    save = TRUE, out = NULL, logger = NULL) {
  # if() stop("not 10X directory")
  if (is.null(sample_index)) sample_index <- basename(d)
  sc <- Read10X(d) %>%
    CreateSeuratObject(
      counts = .,
      project = sample_index,
      min.cells = num_min_cell
    ) %>%
    .reindex()

  if (all(save, !is.null(out))) {
    saveRDS(sc, file = paste0(file.path(out, sample_index), ".RDS"))
    .logger(logger, paste0(": ", sample_index, ".RDS saved"))
  }
  sc
}
.readmtx <- function(
    d, sample_index = NULL, num_min_cell = 3,
    save = TRUE, out = NULL, logger = NULL) {
  # if() stop("not count matirx")
  ## the third line
  if (is.null(sample_index)) sample_index <- basename(d)
  thirdLines <- .readLineN(d, N = 3)
  delim <- SepGuess(thirdLines)

  ## remove duplicated
  df <- read.delim(d, header = TRUE, sep = delim)
  df <- .removeDuplicated(df)
  df <- data.frame(df, row.names = 1)

  sc <- df %>%
    base::as.matrix(x = .) %>%
    as(object = ., "dgCMatrix") %>%
    CreateSeuratObject(
      counts = .,
      project = sample_index,
      min.cells = num_min_cell
    ) %>%
    .reindex()
  if (all(save, !is.null(out))) {
    saveRDS(sc, file = paste0(file.path(out, sample_index), ".RDS"))
    .logger(logger, paste0(": ", sample_index, ".RDS saved"))
  }
  sc
}
.read_h5 <- function(
    d, sample_index = NULL, num_min_cell = 3,
    save = TRUE, out = NULL, logger = NULL) {
  if (is.null(sample_index)) sample_index <- basename(d)
  sc <- Read10X_h5(d) %>%
    CreateSeuratObject(
      counts = .,
      project = sample_index,
      min.cells = num_min_cell
    ) %>%
    .reindex()

  if (all(save, !is.null(out))) {
    saveRDS(sc, file = paste0(file.path(out, sample_index), ".RDS"))
    .logger(logger, paste0(": ", sample_index, ".RDS saved"))
  }
  sc
}

SepGuess <- function(line = NULL) {
  common_sep <- c(",", " ", "\t") # default sep in txt file
  common_sep_tdf <- sapply(
    common_sep, function(x) length(unlist(gregexpr(x, line)))
  )
  if (max(common_sep_tdf) < 2) stop("The frequency of sep Must > 1")
  return(names(which.max(common_sep_tdf)))
}

.readLineN <- function(con = NULL, N = 2) {
  con <- file(con, "r")
  lineN <- readLines(con, n = N)
  return(lineN[N])
  close(con)
}

.removeDuplicated <- function(df = NULL) {
  colnames(df)[1] <- "x"
  df %<>% group_by(x) %>% slice(1)
  # do(head(.,1))
  df
}

.reindex <- function(sc) {
  sampleid <- as.character(sc@meta.data$orig.ident)[1]
  sc <- RenameCells(sc, add.cell.id = sampleid)
  sc
}

runPipe <- function(
    data_combined = NULL, method = c("fastmnn", "harmony"),
    MT_percent = 50,
    nFeature_max = 6000, nFeature_min = 200, resolution = 1,
    logger = NULL, out = "output", doublet = NULL) {
  methods <- c("harmony", "fastmnn")

  if (!any(method %in% methods)) {
    stop("method should be in ", paste(methods, collapse = " "))
  }

  if (!is.null(doublet)) {
    .logger(logger, ": Start removing the doublet cells")
    doub <- read.table(doublet, header = F)
    data_combined <- subset(
      data_combined,
      cells = setdiff(rownames(data_combined@meta.data), doub$V1)
    )
  }

  .logger(logger, ": Start analyzing Mito percentage")
  ## Remove low quality cells, normalization and select top 2,000 variably expressed genes #nolint
  data_combined[["percent.mt"]] <- PercentageFeatureSet(
    data_combined,
    pattern = "^MT-"
  )
  saveRDS(data_combined, "merge_part3.Rds")

  .logger(logger, ": Start preprocessing data.(subseting according to MT, mFeature, Normalization, and FindVariableFeatures)") # nolint
  data_combined %<>% subset(subset = nFeature_RNA > nFeature_min & nFeature_RNA < nFeature_max & percent.mt < MT_percent) %>% # nolint
    NormalizeData(normalization.method = "LogNormalize") %>%
    FindVariableFeatures(selection.method = "vst")

  .logger(logger, ": Scale the data")
  ## Scale the data
  all.genes <- rownames(data_combined)
  data_combined %<>% ScaleData(features = all.genes) %>%
    ScaleData(vars.to.regress = "percent.mt")

  .logger(logger, ": RunPCA for reducing the dimensionality")
  ## RunPCA for reducing the dimensionality
  data_combined %<>% RunPCA(object = ., features = VariableFeatures(object = .))

  str_prefix <- stringi::stri_split(
    as.character(data_combined$orig.ident[1]),
    regex = "_"
  )[[1]][1]

  .logger(logger, ": Start fastmnn Analysis.")
  if ("fastmnn" %in% method) {
    fastmnn_dat <- data_combined
    fastmnn_dat <- # nolint
      RunFastMNN(object.list = SplitObject(fastmnn_dat, split.by = "orig.ident"))

    ## batch effect correction and findclusters

    fastmnn_dat %<>% FindNeighbors(reduction = "mnn", dims = 1:30) %>%
      FindClusters(resolution = resolution)

    fastmnn_dat %<>% RunTSNE(reduction = "mnn", dims = 1:30) %>%
      RunUMAP(reduction = "mnn", dims = 1:30)

    ## save RDS file for subsequent analysis
    saveRDS(
      fastmnn_dat, file.path(out, paste0(str_prefix, "_fastmnn_input.Rds"))
    )
    research_prefix <- paste0(str_prefix, "_fastmnn")
    tsne_umap_save(fastmnn_dat, file.path(out, research_prefix))
  }

  .logger(logger, ": Fastmnn Analysis is Ok -.-")

  ## Calculate the markers for each cluster
  .logger(logger, ": FindMarkers -_-!")
  markers <- FindAllMarkers(fastmnn_dat, only.pos = TRUE, min.pct = 0.25)
  write.table(
    markers, file.path(out, paste0(str_prefix, "_fastmnn_All_Marker.txt")),
    quote = F, sep = "\t"
  )
  .logger(logger, ": FindMarkers is OK -_-!")
  .logger(logger, ": Start harmony Analysis.")
  if ("harmony" %in% method) {
    harmony_dat <- data_combined
    harmony_dat <- harmony_dat %>% RunHarmony("orig.ident")

    harmony_dat %<>% FindNeighbors(reduction = "harmony", dims = 1:30) %>%
      FindClusters(resolution = resolution)

    harmony_dat %<>% RunTSNE(reduction = "harmony", dims = 1:30) %>%
      RunUMAP(reduction = "harmony", dims = 1:30)

    ## save RDS file for subsequent analysis
    saveRDS(
      harmony_dat, file.path(out, paste0(str_prefix, "_harmony_input.Rds"))
    )

    research_prefix <- paste0(str_prefix, "_harmony")
    tsne_umap_save(harmony_dat, file.path(out, research_prefix))
  }
  .logger(logger, ": Harmony Analysis is Ok -.-!")

  ## Calculate the markers for each cluster
  .logger(logger, ": FindMarkers -_-!")
  markers <- FindAllMarkers(harmony_dat, only.pos = TRUE, min.pct = 0.25)
  write.table(
    markers, file.path(out, paste0(str_prefix, "_harmony_All_Marker.txt")),
    quote = F, sep = "\t"
  )
  .logger(logger, ": FindMarkers is OK -_-!")
}

tsne_umap_save <- function(sc, research_prefix) {
  ## umap
  pdf(paste0(research_prefix, "_umap.pdf"), width = 9, height = 9)
  p <- DimPlot(sc, reduction = "umap", label = TRUE, raster = TRUE) + NoLegend()
  print(p)
  dev.off()
  # tsne
  pdf(paste0(research_prefix, "_tsne.pdf"), width = 9, height = 9)
  p <- DimPlot(sc, reduction = "tsne", label = TRUE, raster = TRUE) + NoLegend()
  print(p)
  dev.off()
}

.calltimer <- function() {
  format(Sys.time(), "%Y-%b-%e %H:%M:%S")
}

.logger <- function(logger, text) {
  if (class(logger) == "logger") {
    log4r::info(logger, text)
  } else {
    print(text)
  }
}

#################################################################


#######################  PARAMETERS   ###########################
#################################################################

option_list <- list(
  make_option(c("--path_10x"),
    type = "character",
    default = NULL, help = "the directory of 10x file", metavar = "character"
  ),
  make_option(c("--path_mtx"),
    type = "character",
    default = NULL, help = "the directory of mtx file", metavar = "character"
  ),
  make_option(c("-o", "--out"),
    type = "character",
    default = "preprocess", help = "output directory[default= %default]",
    metavar = "character"
  ),
  make_option(c("--MT_percent"),
    type = "integer",
    default = 50, help = "max percentage of MT", metavar = "number"
  ),
  make_option(c("--min_nFeature"),
    type = "integer",
    default = 200, help = "min nFeature per cell", metavar = "number"
  ),
  make_option(c("--max_nFeature"),
    type = "integer",
    default = 6000, help = "max nFeature per cell", metavar = "number"
  ),
  make_option(c("--resolution"),
    type = "double",
    default = 1, help = "resultion for cluster finder", metavar = "number"
  ),
  make_option(c("--doublet"),
    type = "character",
    default = NULL, help = "the file of doublet", metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

path_10x <- opt[["path_10x"]]
path_mtx <- opt[["path_mtx"]]
percent_MT <- opt[["MT_percent"]]
max_nfeature <- opt[["max_nFeature"]]
min_nfeature <- opt[["min_nFeature"]]
resolution <- opt[["resolution"]]
output <- opt[["out"]]
cell_doublet <- opt[["doublet"]]
##############################################################################

################################# MAIN #####################################
## logger initial


logger <- create.logger()
logfile(logger) <- file.path(output, "log.txt")
level(logger) <- "INFO"


if (!dir.exists(output)) dir.create(output)

res <- NULL
if (!is.null(path_10x)) {
  .logger(logger, ": Reading 10X Data.")
  list_path_10x <- list.dirs(path_10x, full.names = TRUE, recursive = FALSE)
  res <- c(res, parallel::mclapply(
    list_path_10x,
    FUN = function(x, out, logger) .read10x(d = x, out = out, logger = logger),
    mc.cores = mc.cores, out = output, logger = logger
  ))
}

if (!is.null(path_mtx)) {
  .logger(logger, ": Reading Matrix Data")
  list_path_mtx <- list.files(path_mtx, full.names = TRUE, recursive = FALSE)
  if (stringr::str_ends(list_path_mtx[1], ".h5")) {
    read_data <- .read_h5
  } else {
    read_data <- .readmtx
  }
  res <- c(res, parallel::mclapply(
    list_path_mtx,
    FUN = function(x, out, logger) read_data(d = x, out = out, logger = logger),
    out = output, logger = logger, mc.cores = mc.cores
  ))
}

## megre
.logger(logger, ": Multiple samples merge into big Project!")
combined <- Reduce(function(x, y) merge(x, y, merge.data = TRUE), res)

## runpipe
.logger(logger, ": Start RunPipe")
runPipe(
  data_combined = combined, MT_percent = percent_MT,
  nFeature_max = max_nfeature, resolution = resolution,
  nFeature_min = min_nfeature, logger = logger,
  out = output, doublet = cell_doublet
)
