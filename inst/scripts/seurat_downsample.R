#!/usr/bin/env Rscript

###############################################################################
# seurat_downsample.R -- randomly downsample a Seurat object and save it
#
# Draws a uniform random subset of cells (by fraction OR absolute count) with a
# fixed seed and writes a NEW object; the input file is never modified. Reads and
# writes .qs2 (qs2::qs_read/qs_save) or .rds/.RDS (readRDS/saveRDS), inferred
# per-file from the extension.
#
# Default output adds a _ds<tag> and writes to the current directory. A fraction
# tags as a percentage, a count tags in thousands with a trailing "k", so the two
# can never collide:
#   --downsample 0.05  ->  seurat.qs2 -> seurat_ds05.qs2
#   --downsample 0.3   ->  seurat.qs2 -> seurat_ds30.qs2
#   --downsample 30    ->  seurat.qs2 -> seurat_ds0.03k.qs2
#   --downsample 50000 ->  seurat.qs2 -> seurat_ds50k.qs2
#
# Subsetting reconstructs (and validates) any FOV/spatial fields. An object
# written under an older SeuratObject will throw "invalid class FOV object" here;
# pass --update to run UpdateSeuratObject() first (see seurat_update_object.R).
###############################################################################

pacman::p_load(optparse)

option_list <- list(
  make_option("--seurat", type = "character", default = "seurat.qs2",
              help = "Path to the input Seurat object (.qs2 or .rds/.RDS) [default: %default]"),
  make_option("--downsample", type = "numeric", default = 0.05,
              help = "Fraction of cells to retain, in (0, 1] or absolute number of cells > 1 [default: %default]"),
  make_option("--seed", type = "integer", default = 1946,
              help = "RNG seed for sampling [default: %default]"),
  make_option("--outfile", type = "character", default = NULL,
              help = "Output path; format inferred from extension [default: input basename with a _ds<tag>, in the current directory]"),
  make_option("--update", action = "store_true", default = FALSE,
              help = "Run UpdateSeuratObject() before subsetting (for objects written under an older SeuratObject) [default: %default]"),
  make_option("--force", action = "store_true", default = FALSE,
              help = "Overwrite --outfile if it already exists [default: %default]")
)

opts <- parse_args(OptionParser(option_list = option_list))

if (!file.exists(opts$seurat)) {
  stop("--seurat file not found: ", opts$seurat)
}

if (opts$downsample <= 0) {
  stop("--downsample must be > 0")
}

print(opts)



# Resolve the output path (and the _ds<tag>) before the expensive load so a
# clobber is caught immediately.
outfile <- opts$outfile
if (is.null(outfile)) {
  base <- tools::file_path_sans_ext(basename(opts$seurat))
  ext <- tools::file_ext(opts$seurat)
  # Counts always carry a "k" and fractions never do, so a fraction and a count
  # that share digits (0.3 vs 30) cannot resolve to the same filename.
  if (opts$downsample <= 1) {
    tag <- sprintf("ds%02d", round(opts$downsample * 100))
  } else {
    tag <- paste0("ds", format(opts$downsample / 1000, drop0trailing = TRUE,
                               scientific = FALSE, trim = TRUE), "k")
  }
  outfile <- paste0(base, "_", tag, ".", ext)
}

if (file.exists(outfile) && !opts$force) {
  stop("--outfile already exists: ", outfile, " (pass --force to overwrite)")
}


################# RUN #######################
pacman::p_load(sct2, Seurat)

message("\nLoading ", opts$seurat)
seurat <- ReadSeurat(opts$seurat)

if (opts$update) {
  message("Running UpdateSeuratObject()")
  seurat <- UpdateSeuratObject(seurat)
}

# DownsampleObject() clamps a count larger than the object and returns it
# unchanged, which here would write a full-size object under a _ds<tag> naming a
# downsample that never happened. Checked after --update so the count is read off
# a valid object. This is the earliest the object size is known: the _ds<tag> was
# resolved before the load and cannot account for it.
n_cells <- ncol(seurat)
if (opts$downsample > n_cells) {
  stop("--downsample ", opts$downsample, " exceeds the object size ", n_cells)
}

seurat = DownsampleObject(seurat, downsample = opts$downsample, seed = opts$seed)
message("\nDownsampled ", n_cells, " -> ", ncol(seurat), " cells (seed ", opts$seed, ")")
print(seurat)

message("\nSaving to ", outfile)
WriteSeurat(seurat, outfile)

message("\nDone.")


cat("\n\n")
devtools::session_info()
