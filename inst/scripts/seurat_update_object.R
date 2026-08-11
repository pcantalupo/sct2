#!/usr/bin/env Rscript


###############################################################################
# update_seurat_object.R -- load a Seurat object, run UpdateSeuratObject(), save
#
# WHY: objects written under an older SeuratObject (e.g. our HTC env, 5.1.0)
# carry FOV/Centroids objects in the old slot layout. Reading them locally under
# a newer SeuratObject (e.g. 5.4.0) fails on the first subset with:
#   "invalid class FOV object: slots in class definition but not in object:
#    'coords_x_orientation', 'misc'"
# because the @version slot is inherited from an upstream UpdateSeuratObject()
# and does NOT reflect the FOV internals. Re-running UpdateSeuratObject() under
# the newer SeuratObject migrates the FOV/Centroids to the current layout.
#
# Reads/writes .qs2 (qs2::qs_read/qs_save) or .rds/.RDS (readRDS/saveRDS);
# format is inferred per-file from the extension, so --outfile can also convert.
###############################################################################

################# Options ######################
pacman::p_load(optparse)

option_list <- list(
  make_option("--seurat", type = "character", default = NULL,
              help = "Path to the input Seurat object (.qs2 or .rds/.RDS) [required]"),
  make_option("--outfile", type = "character", default = NULL,
              help = "Path to write the updated object; format inferred from extension. [default: overwrite --seurat in place]")
)

opts <- parse_args(OptionParser(option_list = option_list))

if (is.null(opts$seurat)) {
  stop("--seurat is required")
}
if (!file.exists(opts$seurat)) {
  stop("--seurat file not found: ", opts$seurat)
}

outfile <- opts$outfile
if (is.null(outfile)) {
  outfile <- opts$seurat
  message("NOTE: --outfile not given; overwriting --seurat in place:\n  ", outfile)
}
if (!tolower(tools::file_ext(outfile)) %in% c("rds", "qs2")) {
  stop("--outfile must end in .rds or .qs2: ", outfile)
}
if (!dir.exists(dirname(outfile))) {
  stop("--outfile directory does not exist: ", dirname(outfile))
}
if (file.access(dirname(outfile), mode = 2) != 0) {
  stop("--outfile directory is not writable: ", dirname(outfile))
}

print(opts)

######################################################


pacman::p_load(sct2, Seurat)


################# RUN #######################
message("\nInstalled SeuratObject: ", packageVersion("SeuratObject"))

message("\nLoading ", opts$seurat)
seurat <- ReadSeurat(opts$seurat)
message("Object @version before update: ", seurat@version)

message("\nRunning UpdateSeuratObject()")
seurat <- UpdateSeuratObject(seurat)
message("Object @version after update:  ", seurat@version)

# No post-update validation here: UpdateSeuratObject() ends by running
# validObject() on the object and on every image, so a bad migration errors above
# and nothing is written. WriteSeurat() is atomic, so a failed write leaves the
# input intact.

message("\nSaving updated object to ", outfile)
WriteSeurat(seurat, outfile)

message("\nDone.")


cat("\n\n")
devtools::session_info()
