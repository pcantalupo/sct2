#!/usr/bin/env Rscript
pacman::p_load('optparse')

##################### Options ########################
option_list=list(
  make_option("--seuratrds",
              default="",
              type="character",
              help="Seurat object RDS or QS2 [required; default: %default]"),
  make_option("--metadata",
              action="store_true",
              default=FALSE,
              help="Show metadata [optional; default: %default]")
)
opt_parser <- OptionParser(option_list=option_list)
opts <- parse_args(opt_parser)

if (opts$seuratrds == "") {
  message("\nError: --seuratrds parameter is required.")
  print_help(opt_parser)
  quit(status=1)
}
####################################################

pacman::p_load(sct2)

seurat <- ReadSeurat(opts$seuratrds)
SeuratInfo(seurat, metadata = opts$metadata)


