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
# format is inferred per-file from the extension, so --output can also convert.
###############################################################################

################# Options ######################
pacman::p_load(optparse)

option_list <- list(
  make_option("--input", type = "character", default = NULL,
              help = "Path to the input Seurat object (.qs2 or .rds/.RDS) [required]"),
  make_option("--output", type = "character", default = NULL,
              help = "Path to write the updated object; format inferred from extension. [default: overwrite --input in place]"),
  make_option("--check", type = "logical", default = TRUE,
              help = "After updating, subset the first 100 cells to confirm the FOV validates [default: %default]")
)

opts <- parse_args(OptionParser(option_list = option_list))

if (is.null(opts$input)) {
  stop("--input is required")
}
if (!file.exists(opts$input)) {
  stop("--input file not found: ", opts$input)
}

output <- opts$output
if (is.null(output)) {
  output <- opts$input
  message("NOTE: --output not given; overwriting --input in place:\n  ", output)
}

print(opts)

######################################################


pacman::p_load(qs2, Seurat)


############### Helper functions ###############
# Pick reader/writer by file extension; both formats supported on either side.
read_obj <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "qs2") {
    return(qs2::qs_read(path))
  } else if (ext == "rds") {
    return(readRDS(path))
  } else {
    stop("Unsupported input extension '.", ext, "' (expected .qs2 or .rds)")
  }
}

write_obj <- function(seurat, path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "qs2") {
    qs2::qs_save(seurat, path)
  } else if (ext == "rds") {
    saveRDS(seurat, path)
  } else {
    stop("Unsupported output extension '.", ext, "' (expected .qs2 or .rds)")
  }
}


################# RUN #######################
message("\nInstalled SeuratObject: ", packageVersion("SeuratObject"))

message("\nLoading ", opts$input)
seurat <- read_obj(opts$input)
if (!inherits(seurat, "Seurat")) {
  stop("Loaded object is not a Seurat object (class: ", paste(class(seurat), collapse = ", "), ")")
}
message("Object @version before update: ", seurat@version)

message("\nRunning UpdateSeuratObject()")
seurat <- UpdateSeuratObject(seurat)
message("Object @version after update:  ", seurat@version)

# Confirm the migration actually fixed the FOV/Centroids slots: a subset is the
# operation that reconstructs (and validates) the FOV, i.e. the one that was failing.
if (opts$check) {
  message("\nValidation: subsetting first 100 cells to exercise FOV reconstruction")
  n <- min(100, ncol(seurat))
  ok <- tryCatch({
    invisible(subset(seurat, cells = head(Cells(seurat), n)))
    TRUE
  }, error = function(e) {
    message("  VALIDATION FAILED: ", conditionMessage(e))
    FALSE
  })
  if (ok) {
    message("  OK: subset validated, FOV slots are current")
  } else {
    stop("UpdateSeuratObject() did not produce a valid object; NOT saving. ",
         "FOV migration likely incomplete -- regenerate upstream under the current SeuratObject.")
  }
}

message("\nSaving updated object to ", output)
write_obj(seurat, output)

message("\nDone.")


cat("\n\n")
devtools::session_info()
