#!/usr/bin/env Rscript

# Purpose: One UMAP DimPlot per --colorby metadata column (maps to group.by in
#          DimPlot). --colorby takes a comma-separated list so a large object is
#          loaded once and every requested plot is written from that one load.

pacman::p_load(optparse)
pdf(NULL)


##################### Options ########################
option_list <- list(
  make_option("--seurat", default = "", type = "character",
              help = "Seurat object (.rds or .qs2) [required]"),
  make_option("--colorby", default = "", type = "character",
              help = "Comma-separated metadata column(s) that specify the coloring of the cells (group.by); one PNG per column [required]"),
  make_option("--reduction", default = "umap", type = "character",
              help = "Reduction to use [default: %default]"),
  make_option("--label", default = TRUE, type = "logical",
              help = "Label the colorby groups on the plot [default: %default]"),
  make_option("--label_size", default = 4, type = "double",
              help = "Size of the group labels [default: %default]"),
  make_option("--repel", default = FALSE, type = "logical", action = "store_true",
              help = "Repel the group labels [default: %default]"),
  make_option("--height", default = 7, type = "double",
              help = "Plot height in inches [default: %default]"),
  make_option("--width", default = 7, type = "double",
              help = "Plot width in inches [default: %default]"),
  make_option("--outdir", default = ".", type = "character",
              help = "Output directory; PNGs are written to <outdir>/plots [default: %default]")
)
opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

# opts$seurat = "rds/seurat_ds05.qs2"
seuratfile <- opts$seurat
reduction  <- opts$reduction
label      <- opts$label
label_size <- opts$label_size
repel      <- opts$repel
height     <- opts$height
width      <- opts$width
outdir     <- opts$outdir

if (is.null(seuratfile) || seuratfile == "" || !file.exists(seuratfile)) {
  print_help(opt_parser)
  stop("--seurat must be an existing Seurat object (.rds or .qs2)")
}

colorby <- trimws(strsplit(opts$colorby, ",", fixed = TRUE)[[1]])
colorby <- colorby[nzchar(colorby)]
if (length(colorby) == 0) {
  print_help(opt_parser)
  stop("--colorby is required")
}

plotsdir <- file.path(outdir, "plots")
dir.create(plotsdir, recursive = TRUE, showWarnings = FALSE)

message("\nArguments:")
print(opts)
message("\nColor-by columns: ", paste(colorby, collapse = ", "))
message("")
#################################################


pacman::p_load(nvutils, sct2, Seurat, scCustomize, gtools, tidyverse)


message("\nLoading Seurat object from ", seuratfile)
seurat <- ReadSeurat(seuratfile)
print(seurat)

if (!reduction %in% Reductions(seurat)) {
  stop("Reduction not found: '", reduction, "'. Available: ",
       paste(Reductions(seurat), collapse = ", "))
}

# Validate every requested column up front so a typo fails before the first
# (slow) plot rather than partway through the loop.
missing_cols <- setdiff(colorby, colnames(seurat[[]]))
if (length(missing_cols) > 0) {
  stop("Metadata column(s) not found in Seurat object: ", paste(missing_cols, collapse = ", "))
}


################## Plot each --colorby column ###################
outputfiles <- character(0)
for (col in colorby) {
  # Order the groups naturally (0,1,2,...,10 not 0,1,10,2) so the legend and the
  # color assignment are stable, but respect an existing factor's level order
  # rather than clobbering it.
  if (!is.factor(seurat[[]][, col])) {
    col_char <- as.character(seurat[[]][, col])
    seurat[[col]] <- factor(col_char, levels = gtools::mixedsort(unique(col_char)))
  }

  message("\nPlotting ", reduction, " colored-by '", col, "'")
  p <- DimPlot_scCustom(seurat, reduction = reduction, group.by = col,
                        label = label, label.size = label_size, repel = repel) +
    labs(title = paste0(toupper(reduction), " colored by ", col))

  outputfile <- file.path(plotsdir, paste0(toupper(reduction), "_colored_by_", col, ".png"))
  message("Saving to ", outputfile)
  ggsave(outputfile, plot = p, height = height, width = width, bg = "white")
  outputfiles <- c(outputfiles, outputfile)
}


################## Summary ###################
message("\nDone. ", length(outputfiles), " plot(s) saved to ", plotsdir, ":")
for (f in outputfiles) {
  message("  ", basename(f))
}


cat("\n\n")
devtools::session_info()
