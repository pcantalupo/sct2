#' @title ReadSeurat
#' @description Read a Seurat object from an RDS or QS2 file.
#' @param path Path to the file (.rds or .qs2)
#' @return A Seurat object
#' @export
#' @importFrom tools file_ext
#' @examples
#'
#' tmp <- tempfile(fileext = ".rds")
#' WriteSeurat(pbmc_small, tmp)
#' seurat <- ReadSeurat(tmp)
#' unlink(tmp)
#'
ReadSeurat <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "qs2") {
    seurat <- qs2::qs_read(path)
  } else if (ext == "rds") {
    seurat <- readRDS(path)
  } else {
    stop("Unsupported extension '.", ext, "' (expected .rds or .qs2)")
  }
  if (!inherits(seurat, "Seurat")) {
    stop("Loaded object is not a Seurat object (class: ",
         paste(class(seurat), collapse = ", "), ")")
  }
  return(seurat)
}

#' @title WriteSeurat
#' @description Write a Seurat object to an RDS or QS2 file.
#' @param seurat A Seurat object
#' @param path Output path (.rds or .qs2)
#' @return Invisibly returns path
#' @export
#' @importFrom tools file_ext
#' @examples
#'
#' tmp <- tempfile(fileext = ".rds")
#' WriteSeurat(pbmc_small, tmp)
#' unlink(tmp)
#'
WriteSeurat <- function(seurat, path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "qs2") {
    qs2::qs_save(seurat, path)
  } else if (ext == "rds") {
    saveRDS(seurat, path)
  } else {
    stop("Unsupported extension '.", ext, "' (expected .rds or .qs2)")
  }
  invisible(path)
}
