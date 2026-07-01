#!/usr/bin/env Rscript

pacman::p_load(optparse)
pdf(NULL)

# Purpose: UMAP DimPlot colored by celltype, but grouped/labeled by a combined
# <celltype>_<cluster> column so every cluster of a celltype shares that celltype's
# color while staying individually labeled. The combined column is built here from
# the --celltype and --cluster metadata columns.

##################### Options ########################
option_list <- list(
  make_option("--seurat",   default = "",                     type = "character", help = "Seurat object (.rds or .qs2) [required]"),
  make_option("--celltype", default = "singleR_cluster_labels", type = "character", help = "Metadata column with celltype labels (drives colors) [default: %default]"),
  make_option("--cluster",  default = "RNA_snn_res.0.8",      type = "character", help = "Metadata column with cluster IDs (combined with celltype for the label) [default: %default]"),
  make_option("--reduction",  default = "umap",      type = "character", help = "Reduction to use [default: %default]"),
  make_option("--repel",    default = FALSE,                  type = "logical",   action = "store_true", help = "Repel the cluster labels [default: %default]"),
  make_option("--height",   default = 7,                      type = "double",    help = "Plot height in inches [default: %default]"),
  make_option("--width",    default = 7,                      type = "double",    help = "Plot width in inches [default: %default]"),
  make_option("--outputfile", default = "",                   type = "character", help = "Output PNG filename [default: UMAP_colored_by_<celltype>_<cluster>.png]"),
  make_option("--outdir",   default = ".",                    type = "character", help = "Output directory [default: %default]")
)
opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

# opts$seuratfile = "rds/seurat_singleR_prediction.qs2"
seuratfile <- opts$seurat
celltype   <- opts$celltype
cluster    <- opts$cluster
height     <- opts$height
width      <- opts$width
outdir     <- opts$outdir
reduction  <- opts$reduction
repel      <- opts$repel

if (is.null(seuratfile) || seuratfile == "" || !file.exists(seuratfile)) {
  print_help(opt_parser)
  stop("--seurat must be an existing Seurat object (.rds or .qs2)")
}

if (is.null(opts$outputfile) || opts$outputfile == "") {
  opts$outputfile <- paste0("UMAP_colored_by_", celltype, "_", cluster, ".png")
}
outputfile <- opts$outputfile

plotsdir <- file.path(outdir, "plots")
dir.create(plotsdir, recursive = TRUE, showWarnings = FALSE)

message("\nArguments:")
print(opts)
message("")
#################################################

pacman::p_load(nvutils, sct2, Seurat, scCustomize, ggplot2, tidyverse)

message("\nLoading Seurat: ", seuratfile)
seurat <- ReadSeurat(seuratfile)
print(seurat)

md <- seurat[[]]
if (!celltype %in% colnames(md)) {
  stop("celltype column not found in metadata: ", celltype)
}
if (!cluster %in% colnames(md)) {
  stop("cluster column not found in metadata: ", cluster)
}

# Build the combined <celltype>_<cluster> column (e.g. "OB_3", "Macro_5")
celltype_cluster <- "celltype_cluster"
seurat[[celltype_cluster]] <- paste(md[[celltype]], md[[cluster]], sep = "_")
md <- seurat[[]]

# One color per celltype (user's standing discrete-color convention)
labels         <- sort(unique(as.character(md[[celltype]])))
celltypecolors <- set_names(colors_polychrome[seq_along(labels)], labels)

# Expand celltype colors to the celltype_cluster values. The mapping is taken
# directly from the metadata (every celltype_cluster value has exactly one
# celltype), so it is robust to celltype names that contain underscores and to
# any cluster-ID format.
lc_map    <- unique(md[, c(celltype_cluster, celltype)])
lc_colors <- set_names(celltypecolors[as.character(lc_map[[celltype]])],
                       lc_map[[celltype_cluster]])

# DimPlot: colored by celltype, grouped/labeled by celltype_cluster.
# Base DimPlot (not DimPlot_scCustom) because `cols` here is a NAMED vector that
# maps colors to group.by levels by name; scCustomize's colors_use is positional.
message("\nDimPlot colored by celltype, labeled by celltype_cluster")
DimPlot(seurat, group.by = celltype_cluster, cols = lc_colors, reduction = reduction,
        label = TRUE, label.size = 4.5, shuffle = TRUE, repel = repel) +
  labs(title = "UMAP colored by Celltype_Cluster",
       subtitle = paste0("Celltype: ", celltype, " | Clusters: ", cluster)) +
  theme(plot.subtitle = element_text(hjust = 0.5))
ggsave(file.path(plotsdir, outputfile), height = height, width = width, bg = "white")

cat("\n\n")
devtools::session_info()
