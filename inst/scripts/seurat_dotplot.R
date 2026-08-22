#!/usr/bin/env Rscript
pacman::p_load(optparse)   # load other libraries after option parsing so there is less delay on command line.

pdf(NULL)

# Define the options
option_list <- list(
  make_option("--seurat", type="character", default=NULL,
              help="Path to the input Seurat object (.qs2 or .rds/.RDS) [required]"),
  make_option("--markers", type="character", default="rds/markers.rds",
              help="Path to the markers RDS [required]"),
  make_option("--idents", type="character", default="RNA_snn_res.0.8",
              help="Metadata field to use to set the Seurat Idents [default: %default]"),
  make_option("--n_top_genes", type="integer", default=5,
              help="Number of top genes to select [default: %default]"),
  make_option("--title", type="character", default="Top 5 up",
              help="Title of the plot [default: %default]"),
  make_option("--labelsize", type="integer", default=8,
              help="Gene label size [default: %default]"),
  make_option("--rotatelabels", type="logical", default=TRUE,
              help="Rotate identity labels? [default: %default]"),
  make_option("--outputfile", type="character", default="",
              help="Output PNG filename, no directory component [default: dotplot_top<n_top_genes>_<idents>.png]"),
  make_option("--outdir", type="character", default="plots",
              help="Output directory [default: %default]"),
  make_option("--width", type="integer", default=7,
              help="Width of the saved image [default: %default]"),
  make_option("--height", type="integer", default=9,
              help="Height of the saved image [default: %default]")
)

# Parse the options
opt_parser <- OptionParser(option_list=option_list)
opts <- parse_args(opt_parser)

# Check if the seurat file exists
if (is.null(opts$seurat) || opts$seurat == "") {
  print_help(opt_parser)
  stop("Error: The --seurat option is required.\n")
} else if (!file.exists(opts$seurat)) {
  print_help(opt_parser)
  stop("Error: The file ", opts$seurat, " does not exist.\n")
}

# Check if the markers file exists
if (is.null(opts$markers) || opts$markers == "") {
  print_help(opt_parser)
  stop("Error: The --markers option is required.\n")
} else if (!file.exists(opts$markers)) {
  print_help(opt_parser)
  stop("Error: The file ", opts$markers, " does not exist.\n")
}

if (grepl("/", opts$outputfile, fixed = TRUE)) {
  stop("--outputfile must be a filename, not a path; use --outdir for the directory")
}

if (is.null(opts$outputfile) || opts$outputfile == "") {
  opts$outputfile = paste0("dotplot_top", opts$n_top_genes, "_", opts$idents, ".png")
}

dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)

pacman::p_load(sct2, Seurat, ggplot2, dplyr)


# Load data
message("\nLoading Seurat: ", opts$seurat)
seurat <- ReadSeurat(opts$seurat)
markers <- readRDS(opts$markers)

# Set cell type identities
Idents(seurat) <- opts$idents
DefaultAssay(seurat) <- "RNA"

# Select top n genes
top <- markers %>%
  arrange(cluster, -avg_log2FC) %>%
  group_by(cluster) %>%
  top_n(opts$n_top_genes, avg_log2FC)


# Create plot and save as PNG
#   need to use "guides" to change the legend title ('name' param in scale_color_gradientn does not work)
message("\nCreating DotPlot")
gg <- DotPlot(seurat, features = unique(top$gene)) + coord_flip() +
  guides(color = guide_colorbar(title = 'Scaled avg expr')) +
  scale_color_gradientn(colors = c("dodgerblue", "yellow", "indianred")) +
  ggtitle(opts$title)

if (!is.null(opts$labelsize)) {  # I think the default ggplot size is 11
  gg <- gg + theme(axis.text.y = element_text(size = opts$labelsize)) #theme(text = element_text(size = 4))
}
if (!is.null(opts$rotatelabels)) {
  gg <- gg + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}
gg
ggsave(file.path(opts$outdir, opts$outputfile), width = opts$width, height = opts$height, bg = "white")


cat("\n\n")
devtools::session_info()
