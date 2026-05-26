#!/usr/bin/env Rscript
pacman::p_load('optparse')

##################### Options ########################
option_list=list(
  make_option("--seuratrds", default="", type="character", help="Seurat object [required; default: %default]"),
  make_option("--outfile", default="", type="character", help="Output filename to save Seurat metadata [required; default: %default]")
)
opt_parser <- OptionParser(option_list=option_list)
opts <- parse_args(opt_parser)

check_params = function(opt, name) {
  if (opt == "") {
    print_help(opt_parser)
    stop(paste("Argument:", name, " is required"))
  }
}

check_params(opts$seuratrds, "--seuratrds")
check_params(opts$outfile, "--outfile")

#####################################################

pacman::p_load(qs2, sct2)

if (grepl("\\.qs2?$", opts$seuratrds, ignore.case = TRUE)) {
  seurat = qs_read(opts$seuratrds)
} else {
  seurat = readRDS(opts$seuratrds)
}
SaveMetadata(seurat, opts$outfile)

cat("\n\n")
devtools::session_info()

