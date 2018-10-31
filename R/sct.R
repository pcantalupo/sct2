# sct

#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' Seurat2SingleCellExperiment
#'
#' This function creates a Bioconductor SingleCellExperiment from a Seurat object. Raw counts are extracted for the cells used in building the Seurat clusters. The raw counts are normalized by 'scater' package. If you get the following warning message you can ignore it (it comes from scater normalize function): In .local(object, ...) : using library sizes as size factors.
#' @param seurat A Seurat object
#' @param clusterlabels Optional character vector that is same length as the number of clusters. For example, if three clusters (0, 1, 2) then clusterlabels might be c("a","b","c").
#' @keywords Seurat SingleCellExperiment Bioconductor
#' @import Seurat SingleCellExperiment scater
#' @return A SingleCellExperiment object
#' @export
#' @examples
#' Seurat2SingleCellExperiment(seurat, c(1,2,3))   # 'seurat' is created through use of Seurat package

Seurat2SingleCellExperiment = function (seurat, clusterlabels = NULL) {
  require(Seurat)
  require(SingleCellExperiment)
  require(scater)

  cells.to.include = seurat@cell.names  # the cells that were used for cluster determination
  counts = as.data.frame(as.matrix(seurat@raw.data[,cells.to.include]))  # The raw data contains all cells, so need to only select for those that were used for cluster determination. Counts are raw UMIFM counts.

  # annotate each cell with cluster identity and optionally cluster labels
  cluster = as.vector(seurat@ident)
  cell.annotation = data.frame(cluster = cluster, row.names = cells.to.include)
  if (is.character(clusterlabels)) {
    clusternames = sapply(cluster, FUN = function (c, names) { return(names[as.integer(c)+1]) }, clusterlabels, USE.NAMES=F)
    cell.annotation$clusternames = clusternames
  }

  # construct the SingleCellExperiment object
  sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts)), colData = cell.annotation)
  #  str(counts(sce))  # raw count data
  sce = normalize(sce)   # scater function to normalize and log the raw count data
  #  str(exprs(sce))   # the normalized log data
  rowData(sce)$feature_symbol <- rownames(sce)     # use gene names as the feature symbols
  return(sce)
}


#' Self_scmapCluster
#'
#' From a SingleCellExperiment object, this function performs a scmap cluster analysis on itself. Useful for validating clusters.
#' @param sce A SingleCellExperiment object (returned object from Seurat2SingleCellExperiment) with colData populated with 1 or more columns.
#' @keywords scmap
#' @import scmap
#' @return A list of scmapCluster() objects. A self Sankey plot will be generated and opened in a browser for every column in colData()
#' @export
#' @examples
#' names = c("Detox", "DNAReplication" "Quiescent")
#' sce = Seurat2SingleCellExperiment(seurat, names) # 'seurat' is created through use of Seurat package
#' Self_scmapCluster(sce)

Self_scmapCluster = function (sce) {
  require(scmap)
  require(SingleCellExperiment)
  sce <- selectFeatures(sce, suppress_plot = T)

  scr = list()
  for (clustercol in colnames(colData(sce))) {
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

#' TwoSample_scmapCluster
#'
#' This function performs a scmap cluster analysis on two Single Cell Experiment objects. colData columns (in the SCE object) will be compared in serial order
#' @param sce1 A SingleCellExperiment object with colData populated with 1 or more columns (such as the returned object from Seurat2SingleCellExperiment function).
#' @param sce2 Same as 'sce1'. Needs to have the same number of column names in 'colData'.
#' @keywords scmap
#' @import scmap
#' @return A list containing two lists ('1vs2' and '2vs1'). Each list contains a set of scmapCluster() objects. Sankey plots for sce1 vs sce2 and Sankey plots for sce2 vs sce1. Each plot will open in a browser.
#' @export
#' @examples
#' names1 = c("Detox", "DNAReplication" "Quiescent")
#' names2 = c("Cellcycle", "Apoptosis" "Quiescent")
#' sce1 = Seurat2SingleCellExperiment(seurat1, names1) # 'seurat1' is created through use of Seurat package
#' sce2 = Seurat2SingleCellExperiment(seurat2, names2)
#' TwoSample_scmapCluster(sce1, sce2)

TwoSample_scmapCluster = function (sce1, sce2) {
  require(scmap)
  require(SingleCellExperiment)
  sce1 <- selectFeatures(sce1, suppress_plot = T)
  sce2 <- selectFeatures(sce2, suppress_plot = T)

  scr = list()
  # sce1 vs sce2
  for (clustercol in colnames(colData(sce1))) {
    sce2 <- indexCluster(sce2, cluster_col = clustercol)  # adds a Metadata slot called scmap_cluster_index for the 500 feature genes
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

  # sce2 vs sce1
  for (clustercol in colnames(colData(sce2))) {
    sce1 <- indexCluster(sce1, cluster_col = clustercol)  # adds a Metadata slot called scmap_cluster_index for the 500 feature genes
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

