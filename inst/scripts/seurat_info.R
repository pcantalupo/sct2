#!/usr/bin/env Rscript
pacman::p_load('optparse')

##################### Options ########################
option_list=list(
  make_option("--seuratrds",
              default="",
              type="character",
              help="Seurat object RDS or QS2 [required; default: %default]"),
  make_option("--metadata",
              default="",
              type="character",
              help="Save metadata to this TSV file [optional; default: %default]")
)
opt_parser <- OptionParser(option_list=option_list)
opts <- parse_args(opt_parser)

if (opts$seuratrds == "") {
  message("\nError: --seuratrds parameter is required.")
  print_help(opt_parser)
  quit(status=1)
}
####################################################

pacman::p_load(qs2, sct2, Seurat)

# Detect file format and read appropriately
if (grepl("\\.qs2?$", opts$seuratrds, ignore.case = TRUE)) {
  seurat = qs_read(opts$seuratrds)
} else {
  seurat = readRDS(opts$seuratrds)
}
SeuratInfo(seurat)

if (opts$metadata != "") {
  SaveMetadata(seurat, opts$metadata)
}


