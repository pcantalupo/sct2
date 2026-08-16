#' @title ValidateMetadataCols
#' @description Stop with a listing of any \code{cols} not present in a Seurat
#'   object's metadata. Intended for plot scripts to validate user-supplied
#'   column names up front, before the (slow) plotting loop.
#' @param seurat A Seurat object
#' @param cols Character vector of metadata column names to check
#' @return Invisibly returns NULL
#' @export
#' @examples
#'
#' ValidateMetadataCols(pbmc_small, c("orig.ident", "groups"))
#'
ValidateMetadataCols <- function(seurat, cols) {
  missing_cols <- setdiff(cols, colnames(seurat[[]]))
  if (length(missing_cols) > 0) {
    stop("Metadata column(s) not found in Seurat object: ", paste(missing_cols, collapse = ", "))
  }
  invisible(NULL)
}

#' @title ValidateReduction
#' @description Stop, listing the available reductions, if \code{reduction}
#'   is not present in a Seurat object.
#' @param seurat A Seurat object
#' @param reduction Name of the reduction to check
#' @return Invisibly returns NULL
#' @export
#' @importFrom SeuratObject Reductions
#' @examples
#'
#' ValidateReduction(pbmc_small, "pca")
#'
ValidateReduction <- function(seurat, reduction) {
  if (!reduction %in% Reductions(seurat)) {
    stop("Reduction not found: '", reduction, "'. Available: ",
         paste(Reductions(seurat), collapse = ", "))
  }
  invisible(NULL)
}
