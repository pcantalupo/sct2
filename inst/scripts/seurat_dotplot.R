#!/usr/bin/env Rscript
pacman::p_load(optparse)   # load other libraries after option parsing so there is less delay on command line.

pdf(NULL)

# Define the options
option_list <- list(
  make_option(c("-s", "--seurat"), type="character", default=NULL,
              help="Path to the input Seurat object (.qs2 or .rds/.RDS) [required]"),
  make_option(c("-m", "--markers"), type="character", default="rds/markers.rds",
              help="Path to the markers RDS [required]"),
  make_option(c("-i", "--idents"), type="character", default="RNA_snn_res.0.8",
              help="Metadata field to use to set the Seurat Idents [default: %default]"),
  make_option(c("-n", "--n_top_genes"), type="integer", default=5,
              help="Number of top genes to select [default: %default]"),
  make_option(c("-t", "--title"), type="character", default="Top 5 up",
              help="Title of the plot [default: %default]"),
  make_option(c("-z", "--labelsize"), type="integer", default=8,
              help="Gene label size [default: %default]"),
  make_option(c("-r", "--rotatelabels"), type="logical", default=TRUE,
              help="Rotate identity labels? [default: %default]"),
  make_option(c("-o", "--output_path"), type="character", default="dotplot_top5up.png",
              help="Path to save the plot [default: %default]"),
  make_option(c("-w", "--width"), type="integer", default=7,
              help="Width of the saved image [default: %default]"),
  make_option(c("-g", "--height"), type="integer", default=9,
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

# Possible code for dynamically setting the output file name
# if (is.null(opts$output_path)) {
#   markerfilebasename = base::basename(tools::file_path_sans_ext(opts$markers))
#   opts$output_path = paste0("dotplot_top", opts$n_top_genes, markerfilebasename, "_", opts$idents, ".png")
# }
# message("Output file name is: ", opts$output_path); message("\n")

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
ggsave(opts$output_path, width = opts$width, height = opts$height, bg = "white")


cat("\n\n")
devtools::session_info()
