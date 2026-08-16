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
ValidateReduction = function(seurat, reduction) {
  if (!reduction %in% Reductions(seurat)) {
    stop("Reduction not found: '", reduction, "'. Available: ",
         paste(Reductions(seurat), collapse = ", "))
  }
  invisible(NULL)
}
