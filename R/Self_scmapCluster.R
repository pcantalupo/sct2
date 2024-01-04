#' Self_scmapCluster
#'
#' From a SingleCellExperiment object, this function performs a scmap cluster
#' analysis on itself. Useful for validating clusters.
#'
#' @param sce A SingleCellExperiment object
#' @param clustercols Metadata column(s) to use for clustering
#'
#' @keywords scmap
#' @import SingleCellExperiment
#' @importFrom scmap selectFeatures indexCluster scmapCluster getSankey
#' @importFrom SummarizedExperiment colData rowData<-
#' @importFrom S4Vectors metadata
#'
#' @return A list of scmapCluster() objects. A self Sankey plot will be generated
#' for each 'clustercols' and automatically opened in a browser
#' @export
#' @examples
#' \dontrun{
#' Self_scmapCluster(sce)
#' }
Self_scmapCluster = function (sce, clustercols = "RNA_snn_res.0.8") {
  rowData(sce)$feature_symbol <- rownames(sce)   # must have 'feature_symbol' in rowData

  sce <- selectFeatures(sce, suppress_plot = T)

  scr = list()
  for (clustercol in clustercols) {
    message("Running self scmapCluster on ", clustercol)
    sce <- indexCluster(sce, cluster_col = clustercol)  # adds a Metadata slot called scmap_cluster_index for the 500 feature genes

    scr[[clustercol]] <- scmapCluster(
      projection = sce,
      index_list = list(
        self = metadata(sce)$scmap_cluster_index
      )
    )
    plot(getSankey(
      colData(sce)[,clustercol],
      scr[[clustercol]]$scmap_cluster_labs[,'self'],
      plot_height = 400)
    )
  }
  return(scr)
}



