pacman::p_load('optparse')

##################### Options ########################
option_list=list(
  make_option("--seuratrds", default="", type="character", help="Seurat object [required; default: %default]"),
  make_option("--outfile", default="", type="character", help="Output filename to save Seurat metadata [required; default: %default]")
)
opt_parser <- OptionParser(option_list=option_list)
opts <- parse_args(opt_parser)

if (opts$seuratrds == "") {
  message("\nError: --seuratrds is required.")
  print_help(opt_parser)
  quit(status=1)
}

if (opts$outfile == "") {
  message("\nError: --outfile is required.")
  print_help(opt_parser)
  quit(status=1)
}


seurat = readRDS(opts$seuratrds)
md = seurat[[]]
write.table(md, opts$outfile, sep = "\t", quote = FALSE, col.names = NA)



