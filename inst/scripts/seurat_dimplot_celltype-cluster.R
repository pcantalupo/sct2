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
  make_option("--downsample", type = "numeric", default = 1,
              help = "Fraction of cells to retain, in (0, 1] or absolute number of cells > 1 [default: %default]"),
  make_option("--seed", type = "integer", default = 1946,
              help = "RNG seed for sampling [default: %default]"),
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

if (opts$downsample <= 0) {
  stop("--downsample must be > 0")
}

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
n_cells <- ncol(seurat)

# Capturing metadata here b/c any celltype that loses all its cells drops out
# and every celltype after it shifts position in colors_polychrome — the same
# celltype gets a different color at different --downsample values or seeds.
md_full <- seurat[[]]
if (!celltype %in% colnames(md_full)) {
  stop("celltype column not found in metadata: ", celltype)
}
if (!cluster %in% colnames(md_full)) {
  stop("cluster column not found in metadata: ", cluster)
}

# Build the combined <celltype>_<cluster> column (e.g. "OB_3", "Macro_5")
celltype_cluster <- "celltype_cluster"
seurat[[celltype_cluster]] <- paste(md_full[[celltype]], md_full[[cluster]], sep = "_")

# Downsample if requested
if (opts$downsample != 1) {
  message("\nDownsampling Seurat...")
  seurat = DownsampleObject(seurat, downsample = opts$downsample, seed = opts$seed)
  message("\nDownsampled ", n_cells, " -> ", ncol(seurat), " cells (seed ", opts$seed, ")")
  print(seurat)
}

md <- seurat[[]]


# One color per celltype (user's standing discrete-color convention)
labels         <- sort(unique(as.character(md_full[[celltype]])))
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

# Title names the reduction and the actual metadata columns. Subtitle carries the
# downsampling note only, so a full-object plot has no subtitle at all.
title = paste0(toupper(reduction), " colored by ", celltype, " | ", cluster)
n_keep = ncol(seurat)
if (n_keep < n_cells) {
  subtitle = paste0("Cells: ", n_keep, " of ", n_cells)
} else {
  subtitle = NULL
}

DimPlot(seurat, group.by = celltype_cluster, cols = lc_colors,
        reduction = reduction, label = TRUE, label.size = 4.5, shuffle = TRUE,
        repel = repel) +
  labs(title = title, subtitle = subtitle) +
  theme(plot.subtitle = element_text(hjust = 0.5))
ggsave(file.path(plotsdir, outputfile), height = height, width = width, bg = "white")

cat("\n\n")
devtools::session_info()

