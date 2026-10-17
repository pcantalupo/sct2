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
  make_option("--palette", default = "ditto", type = "character",
              help = "Colors for --colorby: ditto (nvutils::colors_ditto, 40 colors), polychrome (nvutils::colors_polychrome, 25 colors) or sccustomize (DimPlot_scCustom default) [default: %default]"),
  make_option("--match", default = "position", type = "character",
              help = "How --palette colors are assigned to --colorby levels: position (nth level in the plot gets the nth color) or name (each factor level keeps one color across runs and subsets) [default: %default]"),
  make_option("--label_size", default = 4, type = "double",
              help = "Size of the group labels [default: %default]"),
  make_option("--repel", default = FALSE, type = "logical", action = "store_true",
              help = "Repel the group labels [default: %default]"),
  make_option("--height", default = 7, type = "double",
              help = "Plot height in inches [default: %default]"),
  make_option("--width", default = 7, type = "double",
              help = "Plot width in inches [default: %default]"),
  make_option("--outdir", default = "plots", type = "character",
              help = "Output directory; one PNG per --colorby column, named <REDUCTION>_colored_by_<colorby>.png [default: %default]")
)
opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

# opts$seurat = "rds/seurat_ds05.qs2"
seuratfile <- opts$seurat
reduction  <- opts$reduction
label      <- opts$label
palette    = opts$palette
match      = opts$match
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

if (!palette %in% c("ditto", "polychrome", "sccustomize")) {
  stop("--palette must be ditto, polychrome or sccustomize")
}
if (!match %in% c("position", "name")) {
  stop("--match must be position or name")
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message("\nArguments:")
print(opts)
message("\nColor-by columns: ", paste(colorby, collapse = ", "))
message("")
#################################################


pacman::p_load(nvutils, sct2, Seurat, scCustomize, gtools, tidyverse)


message("\nLoading Seurat object from ", seuratfile)
seurat <- ReadSeurat(seuratfile)
print(seurat)

ValidateReduction(seurat, reduction)

# Validate every requested column up front so a typo fails before the first
# (slow) plot rather than partway through the loop.
ValidateMetadataCols(seurat, colorby)


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

  # position: the nth level in the plot gets the nth color. name: colors are
  # assigned over the full factor level set and named by level, so a level
  # keeps its color when another level is absent. sccustomize with position
  # passes NULL so DimPlot_scCustom picks its own palette, as before --palette.
  if (match == "name") {
    color_levels = levels(seurat[[]][, col])
  } else {
    color_levels = levels(droplevels(seurat[[]][, col]))
  }
  if (palette == "sccustomize" && match == "position") {
    colors_use = NULL
  } else if (palette == "sccustomize") {
    colors_use = scCustomize_Palette(num_groups = length(color_levels), ggplot_default_colors = FALSE)
  } else if (palette == "polychrome") {
    colors_use = rep(colors_polychrome, length.out = length(color_levels))
  } else {
    colors_use = rep(colors_ditto, length.out = length(color_levels))
  }
  if (match == "name") {
    names(colors_use) = color_levels
  }

  message("\nPlotting ", reduction, " colored-by '", col, "'")
  p <- DimPlot_scCustom(seurat, reduction = reduction, group.by = col,
                        label = label, label.size = label_size, repel = repel,
                        colors_use = colors_use) +
    labs(title = paste0(toupper(reduction), " colored by ", col))

  outputfile <- file.path(outdir, paste0(toupper(reduction), "_colored_by_", col, ".png"))
  message("Saving to ", outputfile)
  ggsave(outputfile, plot = p, height = height, width = width, bg = "white")
  outputfiles <- c(outputfiles, outputfile)
}


################## Summary ###################
message("\nDone. ", length(outputfiles), " plot(s) saved to ", outdir, ":")
for (f in outputfiles) {
  message("  ", basename(f))
}


cat("\n\n")
devtools::session_info()
