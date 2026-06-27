#!/usr/bin/env Rscript

###############################################################################
# seurat_downsample.R -- randomly downsample a Seurat object and save it
#
# Draws a uniform random subset of cells (by fraction OR absolute count) with a
# fixed seed and writes a NEW object; the input file is never modified. Reads and
# writes .qs2 (qs2::qs_read/qs_save) or .rds/.RDS (readRDS/saveRDS), inferred
# per-file from the extension.
#
# Default output adds a _ds<tag> next to the input:
#   --downsample 0.05  ->  seurat.qs2 -> seurat_ds05.qs2
#   --ncells 50000     ->  seurat.qs2 -> seurat_ds50k.qs2
#
# Subsetting reconstructs (and validates) any FOV/spatial fields. An object
# written under an older SeuratObject will throw "invalid class FOV object" here;
# pass --update to run UpdateSeuratObject() first (see seurat_update_object.R).
###############################################################################

pacman::p_load(optparse)

option_list <- list(
  make_option("--seurat", type = "character", default = "seurat.qs2",
              help = "Path to the input Seurat object (.qs2 or .rds/.RDS) [default: %default]"),
  make_option("--downsample", type = "numeric", default = NULL,
              help = "Fraction of cells to retain, in (0, 1]; mutually exclusive with --ncells [default: 0.05 when neither given]"),
  make_option("--ncells", type = "integer", default = NULL,
              help = "Absolute number of cells to retain; mutually exclusive with --downsample [default: %default]"),
  make_option("--seed", type = "integer", default = 1976,
              help = "RNG seed for sampling [default: %default]"),
  make_option("--output", type = "character", default = NULL,
              help = "Output path; format inferred from extension [default: input with a _ds<tag>]"),
  make_option("--update", action = "store_true", default = FALSE,
              help = "Run UpdateSeuratObject() before subsetting (for objects written under an older SeuratObject) [default: %default]"),
  make_option("--force", action = "store_true", default = FALSE,
              help = "Overwrite --output if it already exists [default: %default]")
)

opts <- parse_args(OptionParser(option_list = option_list))

if (!file.exists(opts$seurat)) {
  stop("--seurat file not found: ", opts$seurat)
}

# Resolve sampling mode: fraction (default) or absolute count, never both.
if (!is.null(opts$downsample) && !is.null(opts$ncells)) {
  stop("Specify only one of --downsample or --ncells")
}
if (!is.null(opts$ncells)) {
  mode <- "count"
  if (opts$ncells < 1) {
    stop("--ncells must be >= 1")
  }
} else {
  mode <- "fraction"
  if (is.null(opts$downsample)) {
    frac <- 0.05
  } else {
    frac <- opts$downsample
  }
  if (frac <= 0 || frac > 1) {
    stop("--downsample must be in (0, 1]")
  }
}

print(opts)


# Resolve the output path (and the _ds<tag>) before the expensive load so a
# clobber is caught immediately.
output <- opts$output
if (is.null(output)) {
  base <- tools::file_path_sans_ext(opts$seurat)
  ext <- tools::file_ext(opts$seurat)
  if (mode == "fraction") {
    tag <- sprintf("ds%02d", round(frac * 100))
  } else if (opts$ncells %% 1000 == 0) {
    tag <- paste0("ds", opts$ncells %/% 1000, "k")
  } else {
    tag <- paste0("ds", opts$ncells)
  }
  output <- paste0(base, "_", tag, ".", ext)
}

if (file.exists(output) && !opts$force) {
  stop("--output already exists: ", output, " (pass --force to overwrite)")
}


################# RUN #######################
pacman::p_load(sct2, Seurat)

message("\nLoading ", opts$seurat)
seurat <- ReadSeurat(opts$seurat)
if (!inherits(seurat, "Seurat")) {
  stop("Loaded object is not a Seurat object (class: ", paste(class(seurat), collapse = ", "), ")")
}

if (opts$update) {
  message("Running UpdateSeuratObject()")
  seurat <- UpdateSeuratObject(seurat)
}

n_total <- ncol(seurat)
if (mode == "fraction") {
  n_keep <- floor(n_total * frac)
} else {
  if (opts$ncells > n_total) {
    message("NOTE: --ncells ", opts$ncells, " exceeds object size ", n_total, "; keeping all cells")
  }
  n_keep <- min(opts$ncells, n_total)
}
if (n_keep < 1) {
  stop("Resolved cells-to-keep is ", n_keep, "; nothing to sample")
}

message("\nDownsampling ", n_total, " -> ", n_keep, " cells (seed ", opts$seed, ")")
set.seed(opts$seed)
keep <- sample(Cells(seurat), size = n_keep, replace = FALSE)
seurat <- subset(seurat, cells = keep)
print(seurat)

message("\nSaving to ", output)
WriteSeurat(seurat, output)

message("\nDone.")


cat("\n\n")
devtools::session_info()
