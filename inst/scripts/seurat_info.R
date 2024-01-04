#!/usr/bin/env Rscript
pacman::p_load('optparse')

##################### Options ########################
option_list=list(
  make_option("--seuratrds", default="", type="character", help="Seurat object [required; default: %default]")
)
opt_parser <- OptionParser(option_list=option_list)
opts <- parse_args(opt_parser)

if (opts$seuratrds == "") {
  message("\nError: --seuratrds parameter is required.")
  print_help(opt_parser)
  quit(status=1)
}

pacman::p_load('Seurat', 'sct2')

seurat = readRDS(opts$seuratrds)
SeuratInfo(seurat)

message("\nIdentLabel:")
FindIdentLabel(seurat)

