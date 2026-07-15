#!/usr/bin/env Rscript

# Purpose: Create a split UMAP DimPlot (one panel per --splitby value) colored by
#          a --colorby metadata column (maps to group.by in DimPlot).

pacman::p_load(optparse)
pdf(NULL)


##################### Options ########################
option_list <- list(
  make_option("--seurat", default = "", type = "character",
              help = "Seurat object (.rds or .qs2) [required]"),
  make_option("--splitby", default = "RNA_snn_res.0.8", type = "character",
              help = "Metadata column that specifies individual UMAP plots for each value [default: %default]"),
  make_option("--colorby", default = "orig.ident", type = "character",
              help = "Metadata column that specifies the coloring of the cells (group.by) [default: %default]"),
  make_option("--reduction", default = "umap", type = "character",
              help = "Reduction to use [default: %default]"),
  make_option("--label", default = TRUE, type = "logical",
              help = "Label the colorby groups on the UMAP [default: %default]"),
  make_option("--label_size", default = 4, type = "double",
              help = "Size of the cluster labels [default: %default]"),
  make_option("--repel", default = FALSE, type = "logical", action = "store_true",
              help = "Repel the cluster labels [default: %default]"),
  make_option("--height", default = 7, type = "double",
              help = "Plot height in inches [default: %default]"),
  make_option("--width", default = 7, type = "double",
              help = "Plot width in inches [default: %default]"),
  make_option("--outputfile", default = "", type = "character",
              help = "Output PNG path [default: plots/UMAP_split-<splitby>_color-<colorby>.png]")
)
opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

# opts$seurat = "rds/seurat_ds05.qs2"
seuratfile <- opts$seurat
splitby    <- opts$splitby
colorby    <- opts$colorby
reduction  <- opts$reduction
label      <- opts$label
label_size <- opts$label_size
repel      <- opts$repel
height     <- opts$height
width      <- opts$width

if (is.null(seuratfile) || seuratfile == "" || !file.exists(seuratfile)) {
  print_help(opt_parser)
  stop("--seurat must be an existing Seurat object (.rds or .qs2)")
}

plotsdir <- "plots"
if (is.null(opts$outputfile) || opts$outputfile == "") {
  opts$outputfile <- file.path(plotsdir, paste0("UMAP_split-", splitby, "_color-", colorby, ".png"))
}
outputfile <- opts$outputfile

dir.create(dirname(outputfile), recursive = TRUE, showWarnings = FALSE)

message("\nArguments:")
print(opts)
message("")
#################################################


pacman::p_load(nvutils, sct2, Seurat, scCustomize, gtools, tidyverse)


message("\nLoading Seurat object from ", seuratfile)
seurat <- ReadSeurat(seuratfile)
print(seurat)

# Validate requested metadata columns exist
missing_cols <- setdiff(c(splitby, colorby), colnames(seurat[[]]))
if (length(missing_cols) > 0) {
  stop("Metadata column(s) not found in Seurat object: ", paste(missing_cols, collapse = ", "))
}

# Coerce colorby to a factor so every split panel shares one color scale and a
# single unified legend. A character colorby makes DimPlot_scCustom build a
# per-panel scale from only the values present in that panel, which yields one
# mismatched legend per panel.
if (!is.factor(seurat[[]][, colorby])) {
  colorby_char <- as.character(seurat[[]][, colorby])
  seurat[[colorby]] <- factor(colorby_char, levels = gtools::mixedsort(unique(colorby_char)))
}

# Order the splitby panels naturally (0,1,2,...,10 not 0,1,10,2), but respect an
# existing factor's level order rather than clobbering it.
if (!is.factor(seurat[[]][, splitby])) {
  splitby_char <- as.character(seurat[[]][, splitby])
  seurat[[splitby]] <- factor(splitby_char, levels = gtools::mixedsort(unique(splitby_char)))
}

message("\nPlotting split-by '", splitby, "' colored-by '", colorby, "'")
p <- DimPlot_scCustom(seurat, reduction = reduction, split.by = splitby,
                      group.by = colorby, label = label, label.size = label_size,
                      repel = repel)

message("Saving to ", outputfile)
ggsave(outputfile, plot = p, height = height, width = width, bg = "white")


cat("\n\n")
devtools::session_info()
