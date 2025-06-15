#' @title FixClusterFactorLevels
#' @description
#' Relevel cluster factor levels so they are in numerical order
#' @param seurat A Seurat object
#' @return A Seurat object
#' @export
#' @examples
#'
#' FixClusterFactorLevels(pbmc_small)
#'
FixClusterFactorLevels <- function(seurat) {
  pattern <- "snn_res"
  cluster_resolutions <- grep(pattern, colnames(seurat@meta.data), value = TRUE)

  for (clusters in cluster_resolutions) {
    origfactor = seurat[[]][,clusters]
    newlevels = sort(as.numeric(levels(origfactor)))
    newfactor = factor(seurat[[]][,clusters], levels = newlevels)

    if (!identical(as.character(origfactor), as.character(newfactor))) {
      stop("Original clusters do not equal releveled new clusters")
    } else {
      seurat[[]][,clusters] = newfactor
      message("Releveling ", clusters, " so they are in numerical order")
    }
  }
  return(seurat)
}
