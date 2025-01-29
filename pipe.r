library(tidyverse)
library(Seurat)
library(BisqueRNA)
library(ggplot2)
library(clusterProfiler)
library(Nebulosa)
library(future)
library(reshape2)
library(SingleCellExperiment)
library(dplyr)
library(ggrepel)
library(patchwork)
library(msigdbr)
library(GSVA)
library(RColorBrewer)
library(ggpubr)
library(ROGUE)
library(plyr)
library(viridis)
library(magrittr)
library(data.table)
library(R.utils)
library(grid)
library(cowplot)
library(tidyverse)
library(dittoSeq)
library(harmony)
library(scRepertoire)
library(ggsci)
library(pheatmap)
library(ggpie)
library(sscVis)
library(alakazam)
library(UpSetR)
library(CytoTRACE)
library(ggforce)
library(tidydr)
library(ggplotify)
library(circlize)
library(escape)
library(GSVA)

setwd("/work/xiaxy/work/hcj_xzn/")

data <- readRDS("output/RA_fastmnn_input.Rds")

sub <- data %>%
    subset(idents = c(0:9, 11:16, 18:25)) %>%
    RunUMAP(reduction = "mnn", dims = 1:30) %>%
    RenameIdents(
        "12" = "c01_Monocyte", "18" = "c02_DOCK4+ macrophage",
        "14" = "c03_IEG macrophage", "7" = "c04_S100A12+ macrophage",
        "2" = "c05_SPP1+ macrophage", "8" = "c05_SPP1+ macrophage",
        "19" = "c05_SPP1+ macrophage", "3" = "c06_C1QC+ macrophage",
        "22" = "c07_C1QC+CD16+ macrophage", "13" = "c08_Neutrophil",
        "23" = "c09_DC", "25" = "c10_cDC1", "16" = "c10_cDC1",
        "15" = "c11_cDC2", "21" = "c12_LAMP3+ DC",
        "4" = "c13_Naive T",
        "0" = "c14_CD8+ Tem1", "1" = "c15_CD8+ Tem2",
        "5" = "c16_IEG CD8+ T", "6" = "c17_CXCL13+CD4+ T",
        "9" = "c18_Treg", "11" = "c19_GNLY+ NK",
        "20" = "c20_B cell", "24" = "c21_Fibroblast"
    ) %>%
    mutate(
        patient = str_split(orig.ident, "_", simplify = TRUE) %>%
            apply(1, function(x) paste(x[1], x[2], sep = "_")),
        treatment = ifelse(str_split(orig.ident, "_", simplify = TRUE)[, 3] == 1, "BT", "AT") %>%
            factor(levels = c("BT", "AT")),
        drug = ifelse(patient %in% c("RA_1", "RA_2", "RA_3", "RA_9", "RA_11"), "Adalimumab", "Tofacitinib"),
        anno1 = active.ident,
        anno2 = str_split(as.character(anno1), "_", simplify = TRUE)[, 1]
    )

saveRDS(sub, "input.rds")
