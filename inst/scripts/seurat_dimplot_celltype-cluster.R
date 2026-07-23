#!/usr/bin/env Rscript

pacman::p_load(optparse)
pdf(NULL)

# Purpose: UMAP DimPlot colored by celltype, but grouped/labeled by a combined
# <celltype>_<cluster> column so every cluster of a celltype shares that celltype's
# color while staying individually labeled. The combined column is built here from
# the --celltype and --cluster metadata columns.
#
# --downsample (fraction) or --ncells (absolute count) draws a uniform random
# subset of cells with a fixed --seed; the object itself is not subset.

##################### Options ########################
option_list <- list(
  make_option("--seurat",   default = "",                     type = "character", help = "Seurat object (.rds or .qs2) [required]"),
  make_option("--celltype", default = "singleR_cluster_labels", type = "character", help = "Metadata column with celltype labels (drives colors) [default: %default]"),
  make_option("--cluster",  default = "RNA_snn_res.0.8",      type = "character", help = "Metadata column with cluster IDs (combined with celltype for the label) [default: %default]"),
  make_option("--reduction",  default = "umap",      type = "character", help = "Reduction to use [default: %default]"),
  make_option("--repel",    default = FALSE,                  type = "logical",   action = "store_true", help = "Repel the cluster labels [default: %default]"),
  make_option("--downsample", default = NULL,                 type = "numeric",   help = "Fraction of cells to plot, in (0, 1]; mutually exclusive with --ncells [default: plot all cells]"),
  make_option("--ncells",   default = NULL,                   type = "integer",   help = "Absolute number of cells to plot; mutually exclusive with --downsample [default: plot all cells]"),
  make_option("--seed",     default = 1976,                   type = "integer",   help = "RNG seed for downsampling [default: %default]"),
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
seed       = opts$seed

if (is.null(seuratfile) || seuratfile == "" || !file.exists(seuratfile)) {
  print_help(opt_parser)
  stop("--seurat must be an existing Seurat object (.rds or .qs2)")
}

# Resolve the downsampling mode: fraction, absolute count, or none. Validated
# here so a bad value fails before the expensive load.
if (!is.null(opts$downsample) && !is.null(opts$ncells)) {
  stop("Specify only one of --downsample or --ncells")
}
if (!is.null(opts$ncells)) {
  dsmode = "count"
  if (opts$ncells < 1) {
    stop("--ncells must be >= 1")
  }
} else if (!is.null(opts$downsample)) {
  dsmode = "fraction"
  if (opts$downsample <= 0 || opts$downsample > 1) {
    stop("--downsample must be in (0, 1]")
  }
} else {
  dsmode = "none"
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

# Downsample the cells that get drawn. DimPlot's `cells` argument is used rather
# than subset() so the object (and any FOV/spatial fields) is left untouched and
# the celltype color mapping below still covers every celltype in the object.
n_total = ncol(seurat)
if (dsmode == "fraction") {
  n_keep = floor(n_total * opts$downsample)
} else if (dsmode == "count") {
  if (opts$ncells > n_total) {
    message("NOTE: --ncells ", opts$ncells, " exceeds object size ", n_total, "; plotting all cells")
  }
  n_keep = min(opts$ncells, n_total)
} else {
  n_keep = n_total
}
if (n_keep < 1) {
  stop("Resolved cells-to-plot is ", n_keep, "; nothing to plot")
}
if (n_keep < n_total) {
  message("\nDownsampling ", n_total, " -> ", n_keep, " cells for plotting (seed ", seed, ")")
  set.seed(seed)
  plotcells = sample(Cells(seurat), size = n_keep, replace = FALSE)
} else {
  plotcells = Cells(seurat)
}

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
# Title names the reduction and the actual metadata columns. Subtitle carries the
# downsampling note only, so a full-object plot has no subtitle at all.
title = paste0(toupper(reduction), " colored by ", celltype, " | ", cluster)
if (n_keep < n_total) {
  subtitle = paste0("Cells: ", n_keep, " of ", n_total)
} else {
  subtitle = NULL
}
DimPlot(seurat, cells = plotcells, group.by = celltype_cluster, cols = lc_colors,
        reduction = reduction, label = TRUE, label.size = 4.5, shuffle = TRUE,
        repel = repel) +
  labs(title = title, subtitle = subtitle) +
  theme(plot.subtitle = element_text(hjust = 0.5))
ggsave(file.path(plotsdir, outputfile), height = height, width = width, bg = "white")

cat("\n\n")
devtools::session_info()
