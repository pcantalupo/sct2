#' TwoSample_scmapCluster
#'
#' This function performs a scmap cluster analysis on two Single Cell Experiment
#' objects.
#'
#' @param sce1 A SingleCellExperiment object
#' @param sce2 Same as sce1. Must have same number of column names in 'colData'.
#' @param clustercols Metadata column(s) to use for clustering
#'
#' @keywords scmap
#' @import SingleCellExperiment
#' @importFrom scmap selectFeatures indexCluster scmapCluster getSankey
#' @importFrom SummarizedExperiment colData rowData<-
#' @importFrom S4Vectors metadata
#'
#' @return A list containing two lists ('1vs2' and '2vs1'). Each list contains
#' a set of scmapCluster() objects. Sankey plots for sce1 vs sce2 and Sankey
#' plots for sce2 vs sce1. Each plot will open in a browser.
#' @export
#' @examples
#' \dontrun{
#' TwoSample_scmapCluster(sce1, sce2)
#' }
TwoSample_scmapCluster = function (sce1, sce2, clustercols = "RNA_snn_res.0.8") {
  # must have 'feature_symbol' in rowData
  rowData(sce1)$feature_symbol <- rownames(sce1)
  rowData(sce2)$feature_symbol <- rownames(sce2)

  sce1 <- selectFeatures(sce1, suppress_plot = T)
  sce2 <- selectFeatures(sce2, suppress_plot = T)

  scr = list()
  # sce1 vs sce2
  for (clustercol in clustercols) {
    message("Running sce1 vs sce2 scmapCluster on ", clustercol)

    # adds a Metadata slot called scmap_cluster_index for the 500 feature genes
    sce2 <- indexCluster(sce2, cluster_col = clustercol)
    scr[['1vs2']][[clustercol]] <- scmapCluster(
      projection = sce1,
      index_list = list(
        sce2 = metadata(sce2)$scmap_cluster_index
      )
    )
    plot(getSankey(
      colData(sce1)[,clustercol],
      scr[['1vs2']][[clustercol]]$scmap_cluster_labs[,'sce2'],
      plot_height = 400)
    )
  }

  cat("\n")

  # sce2 vs sce1
  for (clustercol in clustercols) {
    message("Running sce2 vs sce1 scmapCluster on ", clustercol)

    # adds a Metadata slot called scmap_cluster_index for the 500 feature genes
    sce1 <- indexCluster(sce1, cluster_col = clustercol)
    scr[['2vs1']][[clustercol]] <- scmapCluster(
      projection = sce2,
      index_list = list(
        sce1 = metadata(sce1)$scmap_cluster_index
      )
    )
    plot(getSankey(
      colData(sce2)[,clustercol],
      scr[['2vs1']][[clustercol]]$scmap_cluster_labs[,'sce1'],
      plot_height = 400)
    )
  }
  return(scr)
}


