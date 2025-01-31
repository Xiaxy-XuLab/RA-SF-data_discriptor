# Load necessary libraries
if (T) {
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
    library(scales)
    library(doParallel)
    library(readr)
    library(future)
    library(org.Hs.eg.db)
}

# Set working directory
setwd("/work/xiaxy/work/hcj_xzn/")

# Load data
data <- readRDS("input.rds")
data <- readRDS("input_del_neu.rds")

# Define color palettes
type_color <- c(
    "#8CC8D6", "#1078a9", "#ECB27F", "#dd626f", "#659E47",
    "#C9B3CF", "#c373ac", "#9EC693", "#E2833B", "#815AA8",
    "#CB3F36", "#A3A1A6", "#D9CE99", "#B294C7", "#EBD57C",
    "#E93A3B", "#C48244", "#D7A7B4", "#506397", "#45676D",
    "#A4A399", "#3ea2e5"
)

group_color <- c("#A46D32", "#32A19C")

ty_color <- c(
    "#8CC8D6", "#dd626f", "#ECB27F", "#1078a9", "#659E47",
    "#815AA8", "#506397", "#EBD57C", "#C9B3CF", "#c373ac",
    "#9EC693", "#E2833B", "#CB3F36", "#A3A1A6", "#D9CE99",
    "#B294C7", "#E93A3B", "#C48244", "#D7A7B4", "#506397",
    "#45676D", "#A4A399", "#3ea2e5"
)

# Plot UMAP for anno3
pdf("Figure1/p1_umap_t.pdf", width = 30, height = 30)
DimPlot(data,
    label = FALSE, group.by = "anno3",
    cols = type_color,
    pt.size = 0.1, raster = FALSE
) +
    NoAxes() + NoLegend() +
    labs(title = "")
dev.off()

# Plot UMAP for treatment
pdf("Figure1/p2_umap_t.pdf", width = 30, height = 30)
DimPlot(data,
    label = FALSE, group.by = "treatment",
    cols = group_color,
    pt.size = 0.1, raster = FALSE
) +
    NoAxes() +
    labs(title = "")
dev.off()

# Calculate and plot fraction barplot
mat <- data@meta.data %>%
    select(treatment, anno3) %>%
    table() %>%
    melt() %>%
    mutate(value = value / as.vector(table(data$treatment)[treatment])) %>%
    mutate(treatment = factor(treatment, levels = c("BT", "AT")))

p <- ggplot(mat, aes(x = treatment, y = value, fill = anno3)) +
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values = ty_color[c(1:5, 7)]) +
    theme_classic() +
    labs(x = "", y = "Fraction") +
    theme(
        legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 12)
    ) +
    scale_y_continuous(expand = c(0, 0.01), labels = percent)
ggsave("Figure1/e_percent.pdf", p, width = 3, height = 2.7)

# Process data and create fraction barplot
mat <- data@meta.data %>%
    select(orig.ident, anno3) %>%
    table() %>%
    melt() %>%
    mutate(value = value / as.vector(table(data$orig.ident)[orig.ident])) %>%
    mutate(
        patient = str_split_fixed(as.character(orig.ident), "_", 3) %>%
            as_tibble() %>%
            unite("patient", V1, V2, sep = "_") %>%
            pull(patient),
        value = ifelse(str_split_fixed(as.character(orig.ident), "_", 3)[, 3] == "2", value, -1 * value) # nolint
    ) %>%
    mutate(
        patient = factor(patient, levels = rev(unique(patient)[c(1, 4:9, 2, 3)])), # nolint
        anno3 = factor(anno3, levels = rev(levels(anno3)))
    )

# Plot fraction barplot
p <- mat %>%
    ggplot(aes(x = value, y = patient, fill = anno3)) +
    geom_bar(position = "fill", stat = "identity", width = 0.85) +
    scale_fill_manual(values = rev(ty_color[c(1:5, 7)])) +
    theme_bw() +
    guides(fill = guide_legend(title = NULL)) +
    labs(y = "", x = "Fraction") +
    theme(
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_text(color = "black", size = 10),
        axis.text.x = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 12)
    ) +
    scale_x_continuous(
        expand = c(0, 0.01),
        labels = percent, position = "top"
    ) +
    guides(fill = guide_legend(reverse = TRUE))

# Save the plot
ggsave("Figure1/d_percent.pdf", p, width = 4.3, height = 3)

# Define genes for violin plot
gene <- c("IL6", "TNF", "IL1A", "IL1B", "IL12A", "IL12B", "IL18", "JAK1", "JAK2", "JAK3") # nolint

# Plot violin plot
pdf("Figure1/vln.pdf", width = 4.5, height = 6)
VlnPlot(
    data,
    features = gene,
    stack = TRUE,
    group.by = "anno3",
    sort = FALSE,
    flip = TRUE,
    cols = ty_color[c(1:5, 7)],
    pt.size = 0,
    split.by = "anno3"
) +
    theme(
        legend.position = "none",
        axis.title.x = element_blank()
    )
dev.off()

# Set up parallel processing
plan("multicore", workers = 60)
options(future.globals.maxSize = 100000 * 1024^2)

# Subset data for Macrophage and T cell analysis
macro_tof <- data %>%
    subset(anno3 == "Macrophage" & drug == "Tofacitinib") %>%
    {
        .@active.ident <- factor(.$treatment)
        .
    }

tcell_ada <- data %>%
    subset(anno3 == "T cell" & drug == "Adalimumab") %>%
    {
        .@active.ident <- factor(.$treatment)
        .
    }

# Find markers for Macrophage and T cell subsets
marker_macro_tof <- FindAllMarkers(macro_tof, only.pos = TRUE, min.pct = 0, logfc.threshold = 0) # nolint
marker_tcell_ada <- FindAllMarkers(tcell_ada, only.pos = TRUE, min.pct = 0, logfc.threshold = 0) # nolint

# Save marker results
write.table(marker_macro_tof, "input/macrophage_deg_BT-AT_tof.txt", quote = FALSE, sep = "\t") # nolint
write.table(marker_tcell_ada, "input/tcell_deg_BT-AT_ada.txt", quote = FALSE, sep = "\t") # nolint

# Define helper functions
gene_subset <- function(x, y) {
    x %>%
        filter(cluster == y & p_val_adj < 10^-50) %>%
        pull(gene)
}

pos_calculate <- function(x) {
    100 * sum(x > 0) / length(x)
}

vol_graph <- function(x, y, z, k, m) {
    x <- x %>%
        mutate(
            avg_log2FC = ifelse(cluster == y, -1 * avg_log2FC, avg_log2FC),
            p_val_adj = ifelse(p_val_adj == 0, min(p_val_adj[p_val_adj != 0]), p_val_adj), # nolint
            label = ifelse(gene %in% z, gene, NA),
            group = case_when(
                p_val_adj >= 10^-10 ~ "N.S.",
                gene %in% k & avg_log2FC > 0 ~ "Up",
                gene %in% k & avg_log2FC < 0 ~ "Down",
                TRUE ~ "Sig"
            ),
            group = factor(group, levels = c("N.S.", "Sig", "Down", "Up"))
        ) %>%
        arrange(group)

    ggplot(x, aes(x = avg_log2FC, y = -log10(p_val_adj), color = group)) +
        geom_point(size = 1) +
        scale_color_manual(values = c("lightgrey", "gray", "#77A270", "#E3B66B")) + # nolint
        scale_y_continuous(limits = c(0, 310)) +
        geom_text_repel(aes(label = label), color = "#3750A1", max.overlaps = 30, size = 2.5) + # nolint
        theme_classic() +
        labs(title = m, y = "", x = "") +
        theme(
            plot.title = element_text(hjust = 0.5, size = 12),
            legend.position = "none"
        )
}

plotenrich <- function(y) {
    plotEnrichment(fgsea_sets[[y]], ranks) +
        labs(title = y, x = "Rank", y = "Enrichment score") +
        theme(
            plot.title = element_text(hjust = 0.5, color = "black", size = 8),
            axis.text = element_text(color = "black")
        )
}

# Read marker files
marker_a <- read.table("input/macrophage_deg_BT-AT.txt", header = TRUE, sep = "\t") # nolint
marker_b <- read.table("input/macrophage_deg_BT-AT_ada.txt", header = TRUE, sep = "\t") # nolint
marker_c <- read.table("input/macrophage_deg_BT-AT_tof.txt", header = TRUE, sep = "\t") # nolint

# Find intersecting genes
inter_neg <- Reduce(intersect, list(
    gene_subset(marker_a, "AT"),
    gene_subset(marker_b, "AT"),
    gene_subset(marker_c, "AT")
))

inter_pos <- Reduce(intersect, list(
    gene_subset(marker_a, "BT"),
    gene_subset(marker_b, "BT"),
    gene_subset(marker_c, "BT")
))

# Define gene labels
label <- c(
    "CXCL8", "SPP1", "CXCL3", "CXCL2", "SLC2A3", "CCL7",
    "VEGFA", "STAT1", "SLAMF9", "JUN", "CCL2", "NFKBIA",
    "GBP5", "CXCL5", "IL1B"
)

# Generate and save volcano plots
final <- vol_graph(marker_a, "AT", label, c(inter_pos, inter_neg), "") +
    vol_graph(marker_b, "AT", label, c(inter_pos, inter_neg), "Adalimumab") +
    vol_graph(marker_c, "AT", label, c(inter_pos, inter_neg), "Tofacitinib") +
    plot_layout(nrow = 1)

ggsave("Figure1/h_vol.pdf", final, width = 12, height = 4)

# Read marker files
marker_files <- c(
    "input/tcell_deg_BT-AT.txt",
    "input/tcell_deg_BT-AT_ada.txt",
    "input/tcell_deg_BT-AT_tof.txt"
)

markers <- marker_files %>%
    map(~ read.table(.x, header = TRUE, sep = "\t"))

# Extract gene subsets for AT and BT clusters
extract_gene_subsets <- function(marker, cluster) {
    marker %>%
        filter(cluster == !!cluster & p_val_adj < 10^-50) %>%
        pull(gene)
}

inter_neg <- markers %>%
    map(~ extract_gene_subsets(.x, "AT")) %>%
    reduce(intersect)

inter_pos <- markers %>%
    map(~ extract_gene_subsets(.x, "BT")) %>%
    reduce(intersect)

# Define gene labels
label <- c(
    "JAK3", "GZMB", "SPP1", "LAG3", "STAT3", "IFI16", "DUSP4",
    "HAVCR2", "NFKBIZ", "GBP1", "IRF3", "CCR5", "MT1E", "MT1X",
    "S100A6", "ANXA1"
)

# Generate volcano plots
volcano_plots <- list(
    vol_graph(markers[[1]], "AT", label, c(inter_pos, inter_neg), ""),
    vol_graph(markers[[2]], "AT", label, c(inter_pos, inter_neg), "Adalimumab"),
    vol_graph(markers[[3]], "AT", label, c(inter_pos, inter_neg), "Tofacitinib")
)

# Combine plots and save
final_plot <- wrap_plots(volcano_plots, nrow = 1)
ggsave("Figure1/h_vol_tcell.pdf", final_plot, width = 12, height = 4)

#############################################################
# Define helper functions
extract_gene_subsets <- function(marker, cluster) {
    marker %>%
        filter(cluster == !!cluster & p_val_adj < 10^-50) %>%
        pull(gene)
}

plot_volcano <- function(marker, cluster, label, inter_genes, title) {
    vol_graph(marker, cluster, label, inter_genes, title)
}

plot_go_enrichment <- function(go_result, select_terms, output_file, width, height) { # nolint
    go_result %>%
        filter(Description %in% select_terms) %>%
        mutate(
            fold = (as.numeric(str_split_fixed(GeneRatio, "/", 2)[, 1]) / as.numeric(str_split_fixed(GeneRatio, "/", 2)[, 2])) / # nolint
                (as.numeric(str_split_fixed(BgRatio, "/", 2)[, 1]) / as.numeric(str_split_fixed(BgRatio, "/", 2)[, 2])), # nolint
            Description = factor(Description, levels = Description[order(fold)])
        ) %>%
        ggplot(aes(x = fold, y = Description)) +
        geom_point(aes(color = p.adjust, size = Count)) +
        scale_color_gradient(low = "red", high = "blue") +
        xlab("Fold Enrichment") +
        theme_bw() +
        theme(
            axis.title.y = element_blank(),
            axis.text = element_text(color = "black")
        )
    ggsave(output_file, width = width, height = height)
}

# Read marker files
marker_files <- c(
    "input/macrophage_deg_BT-AT.txt",
    "input/macrophage_deg_BT-AT_ada.txt",
    "input/macrophage_deg_BT-AT_tof.txt"
)

markers <- marker_files %>%
    map(~ read.table(.x, header = TRUE, sep = "\t"))

# Extract gene subsets for AT and BT clusters
inter_neg <- markers %>%
    map(~ extract_gene_subsets(.x, "AT")) %>%
    reduce(intersect)

inter_pos <- markers %>%
    map(~ extract_gene_subsets(.x, "BT")) %>%
    reduce(intersect)

# Define gene labels
label <- c(
    "CXCL8", "SPP1", "CXCL3", "CXCL2", "SLC2A3", "CCL7",
    "VEGFA", "STAT1", "SLAMF9", "JUN", "CCL2", "NFKBIA",
    "GBP5", "CXCL5", "IL1B"
)

# Generate volcano plots
volcano_plots <- list(
    plot_volcano(markers[[1]], "AT", label, c(inter_pos, inter_neg), ""),
    plot_volcano(markers[[2]], "AT", label, c(inter_pos, inter_neg), "Adalimumab"), # nolint
    plot_volcano(markers[[3]], "AT", label, c(inter_pos, inter_neg), "Tofacitinib") # nolint
)

# Combine plots and save
final_plot <- wrap_plots(volcano_plots, nrow = 1)
ggsave("Figure1/h_vol.pdf", final_plot, width = 12, height = 4)

# Perform GO enrichment analysis
encirh_go_bp <- enrichGO(
    gene = inter_pos,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
)

# Define selected GO terms
select_terms <- c(
    "response to hypoxia", "regulation of angiogenesis",
    "response to lipopolysaccharide", "cytokine-mediated signaling pathway",
    "response to interferon-gamma", "cellular response to interleukin-1",
    "response to tumor necrosis factor"
)

# Plot and save GO enrichment results
plot_go_enrichment(encirh_go_bp@result, select_terms, "Figure1/macro_go.pdf", 6, 3.5) # nolint

# Repeat for T cell data
marker_files_tcell <- c(
    "input/tcell_deg_BT-AT.txt",
    "input/tcell_deg_BT-AT_ada.txt",
    "input/tcell_deg_BT-AT_tof.txt"
)

markers_tcell <- marker_files_tcell %>%
    map(~ read.table(.x, header = TRUE, sep = "\t"))

inter_neg_tcell <- markers_tcell %>%
    map(~ extract_gene_subsets(.x, "AT")) %>%
    reduce(intersect)

inter_pos_tcell <- markers_tcell %>%
    map(~ extract_gene_subsets(.x, "BT")) %>%
    reduce(intersect)

# Generate volcano plots for T cell data
volcano_plots_tcell <- list(
    plot_volcano(markers_tcell[[1]], "AT", label, c(inter_pos_tcell, inter_neg_tcell), ""), # nolint
    plot_volcano(markers_tcell[[2]], "AT", label, c(inter_pos_tcell, inter_neg_tcell), "Adalimumab"), # nolint
    plot_volcano(markers_tcell[[3]], "AT", label, c(inter_pos_tcell, inter_neg_tcell), "Tofacitinib") # nolint
)

# Combine plots and save
final_plot_tcell <- wrap_plots(volcano_plots_tcell, nrow = 1)
ggsave("Figure1/h_vol_tcell.pdf", final_plot_tcell, width = 12, height = 4)

# Perform GO enrichment analysis for T cell data
encirh_go_bp_tcell <- enrichGO(
    gene = inter_pos_tcell,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
)

# Define selected GO terms for T cell data
select_terms_tcell <- c(
    "regulation of innate immune response", "T cell activation",
    "tumor necrosis factor-mediated signaling pathway",
    "T cell receptor signaling pathway",
    "response to type I interferon",
    "response to oxygen levels",
    "regulation of I-kappaB kinase/NF-kappaB signaling",
    "positive regulation of cytokine production",
    "T cell differentiation"
)

# Plot and save GO enrichment results for T cell data
plot_go_enrichment(encirh_go_bp_tcell@result, select_terms_tcell, "Figure1/tcell_go.pdf", 6, 4) # nolint
#############################################################
