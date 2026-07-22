#!/usr/bin/env Rscript


###############################################################################
# seurat_strip_scaledata.R -- drop scale.data layers from every assay and save
#
# WHY: scale.data is a dense matrix (features x cells) and usually dominates the
# on-disk size of a saved object, yet it is trivially regenerated with
# ScaleData(). Stripping it makes objects far cheaper to move and load.
#
# Handles both assay classes: v5 assays (Assay5) can hold several scale.data
# layers (scale.data, scale.data.1, ...), all of which are removed; v3 assays
# (Assay) hold a single @scale.data slot, which is reset to an empty matrix.
#
# Overwrites --seurat in place by default; pass --outfile to write elsewhere.
# Reads/writes .qs2 (qs2::qs_read/qs_save) or .rds/.RDS (readRDS/saveRDS);
# format is inferred per-file from the extension, so --outfile can also convert.
###############################################################################

################# Options ######################
pacman::p_load(optparse)

option_list = list(
  make_option("--seurat", type = "character", default = "seurat.qs2",
              help = "Path to the input Seurat object (.qs2 or .rds/.RDS) [default: %default]"),
  make_option("--outfile", type = "character", default = NULL,
              help = "Path to write the stripped object; format inferred from extension [default: overwrite --seurat in place]"),
  make_option("--force", action = "store_true", default = FALSE,
              help = "Overwrite --outfile if it already exists [default: %default]")
)

opts = parse_args(OptionParser(option_list = option_list))

if (!file.exists(opts$seurat)) {
  stop("--seurat file not found: ", opts$seurat)
}

inplace = is.null(opts$outfile)
output = opts$outfile
if (inplace) {
  output = opts$seurat
  message("NOTE: --outfile not given; overwriting --seurat in place:\n  ", output)
}

# Only guard a distinct --outfile; an in-place rewrite is the documented default.
if (!inplace && file.exists(output) && !opts$force) {
  stop("--outfile already exists: ", output, " (pass --force to overwrite)")
}

print(opts)

######################################################


pacman::p_load(sct2, Seurat)


################# RUN #######################
message("\nLoading ", opts$seurat)
seurat = ReadSeurat(opts$seurat)
message("Object size before: ", format(object.size(seurat), units = "auto"))

# Counted separately: a v5 assay can hold several scale.data layers, so layers
# and assays are not interchangeable units.
n_layers = 0
n_assays = 0

for (assay in Assays(seurat)) {
  assay_obj = seurat[[assay]]

  if (inherits(assay_obj, "Assay5")) {
    # v5: scale.data is one or more named layers; drop each by assigning NULL.
    layers = grep("^scale\\.data", Layers(assay_obj), value = TRUE)
    for (layer in layers) {
      message("Assay ", assay, ": dropping layer '", layer, "'")
      LayerData(seurat[[assay]], layer = layer) = NULL
      n_layers = n_layers + 1
    }
    if (length(layers) > 0) {
      n_assays = n_assays + 1
    }
  } else {
    # v3: scale.data is a fixed slot, so reset it to the empty default rather
    # than deleting it.
    if (nrow(assay_obj@scale.data) > 0) {
      message("Assay ", assay, ": clearing @scale.data slot (",
              nrow(assay_obj@scale.data), "x", ncol(assay_obj@scale.data), ")")
      assay_obj@scale.data = new(Class = "matrix")
      seurat[[assay]] = assay_obj
      n_layers = n_layers + 1
      n_assays = n_assays + 1
    }
  }
}

if (n_assays == 0) {
  message("\nNo scale.data found in any assay; object unchanged")
} else {
  message("\nDropped ", n_layers, " scale.data layer(s) from ", n_assays, " assay(s)")
}
message("Object size after:  ", format(object.size(seurat), units = "auto"))

print(seurat)

message("\nSaving to ", output)
WriteSeurat(seurat, output)

message("\nDone.")


cat("\n\n")
devtools::session_info()
