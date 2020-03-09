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
#' \dontrun{
#' names = c("Detox", "DNAReplication", "Quiescent")
#' sce = Seurat2SingleCellExperiment(seurat, names) # 'seurat' is created through use of Seurat package
#' Self_scmapCluster(sce)
#' }

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
#' names1 = c("Detox", "DNAReplication", "Quiescent")
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


#' addBKVFracExpression_to_Seurat
#'
#' This function adds seven columns to slot @meta.data of a Seurat object. They are 1) Total_mRNA, 2) counts_VP1, 3) counts_LTAg, 4) counts_BKV, 5) frac_VP1, 6) frac_LTAg, and 7) frac_BKV. 'Total_mRNA' is the sum of the raw UMI counts in each cell. The 'counts' are sum of raw counts for VP1 or LTAg in each cell. 'Frac' is the fraction of VP1 (or LTAg, BKV) counts in each cell using 'Total_mRNA'.
#' @param seurat A Seurat object
#' @keywords Seurat
#' @import Seurat
#' @return A Seurat object with 7 new columns in @meta.data
#' @export
#' @examples
#' seurat = addBKVFracExpression_to_Seurat(seurat)

addBKVFracExpression_to_Seurat = function(seurat) {
  require(Seurat)
  rawdata.filt = seurat@raw.data[,seurat@cell.names]  # get raw data for the filtered cells
  bkv.data = data.frame("Total_mRNA" = colSums(as.matrix(rawdata.filt)))

  bkv.genes=c("VP1", "LTAg")  
  # create data.frame of all zeros for default gene expression
  tmp = data.frame(matrix(data = 0, ncol = length(bkv.genes), nrow = nrow(bkv.data)))
  colnames(tmp) = bkv.genes
  # fill data.frame with gene expression values
  for (gene in bkv.genes) {
    if (gene %in% rownames(rawdata.filt)) {  # check if bkv gene is found
      tmp[,gene] = rawdata.filt[gene, ]
    }
  }

  bkv.data = cbind(bkv.data, tmp)
  bkv.data$BKV = rowSums(bkv.data[,bkv.genes])
  colnames(bkv.data) = c("Total_mRNA", "counts_VP1", "counts_LTAg", "counts_BKV")
  bkv.data$frac_VP1 = bkv.data$counts_VP1/bkv.data$Total_mRNA
  bkv.data$frac_LTAg = bkv.data$counts_LTAg/bkv.data$Total_mRNA
  bkv.data$frac_BKV = bkv.data$counts_BKV/bkv.data$Total_mRNA
  seurat = AddMetaData(seurat, metadata = bkv.data)
  return(seurat)
}


#' @export
# seurat A Seurat object
# name The meta.data field that needs releveling
# groups i.e. c("0%", "1%", "10%", "100%")
# bkvname either of c("BKVm5", "LT")
relevel_groupgenotype_in_Seurat = function (seurat, name, groups, bkvname) {
  if (missing(seurat) || missing(name) || missing(groups) || missing(bkvname)) {
    stop("Must supply seurat, name, groups, and bkvname to 'myrelevel'")
  }
  v = seurat@meta.data[,name]
  mylevels = c(paste0(groups, "_Mock"), paste0(groups, paste0("_", bkvname)))
  mymatch = match(mylevels, levels(v), nomatch = 0)
  seurat@meta.data[,name] = factor(v, levels=levels(v)[mymatch])
  seurat
}

#' addBKVExprGroups_to_Seurat
#'
#' This function adds 4 columns to slot @meta.data of a Seurat object. They are 1) VP1_group, 2) VP1_group_genotype, 3) LTAg_group, and 4) LTAg_group_genotype.
#' @param seurat A Seurat object
#' @keywords Seurat
#' @import Seurat
#' @return A Seurat object with 4 new columns in @meta.data
#' @export
#' @examples
#' seurat = addBKVExprGroups_to_Seurat(seurat)
addBKVExprGroups_to_Seurat = function(seurat) {
  require(Seurat)
  require(dplyr)

  # Add VP1 group info to Seurat
  message("Adding VP1 group metadata")
  VP1_groups = c("0%", "<=1%", "<=10%", "<=100%")
  c = as.character(cut(seurat@meta.data$frac_VP1, breaks=c(0, 0.01, 0.1, 1)))
  c[is.na(c)] = "0%"
  c = factor(recode(c, "(0,0.01]" = "<=1%", "(0.01,0.1]" = "<=10%", "(0.1,1]" = "<=100%"), levels = VP1_groups)
  c2 = factor(paste0(c, "_", seurat@meta.data$genotype))#, levels=paste0(VP1_groups, "_BKVm5"))
  seurat@meta.data$VP1_group = c
  seurat@meta.data$VP1_group_genotype = c2
  seurat = relevel_groupgenotype_in_Seurat(seurat, name = "VP1_group_genotype", groups = VP1_groups, bkvname = "BKVm5")

  # table(seurat@meta.data$VP1_group, useNA = "always")
  # table(seurat@meta.data$VP1_group_genotype, useNA = "always")
  
  # Add LTAg group info to Seurat (max 0%, 0.05% 1:2000, 0.2% 1:500, 10% 1:10)
  message("Adding LTAg group metadata")
  LTAg_groups = c("0%", "<=0.05%", "<=0.2%", "<=10%")
  c = as.character(cut(seurat@meta.data$frac_LTAg, breaks=c(0, 0.0005, 0.002, 0.1)))
  c[is.na(c)] = "0%"
  c = factor(recode(c, "(0,0.0005]" = "<=0.05%", "(0.0005,0.002]" = "<=0.2%", "(0.002,0.1]" = "<=10%"), levels = LTAg_groups)
  c2 = factor(paste0(c, "_", sub("BKVm5", "LT", seurat@meta.data$genotype)))#, levels=paste0(LTAg_groups, "_LT"))
  seurat@meta.data$LTAg_group = c
  seurat@meta.data$LTAg_group_genotype = c2
  seurat = relevel_groupgenotype_in_Seurat(seurat, name = "LTAg_group_genotype", groups = LTAg_groups, bkvname = "LT")
  
  # table(seurat@meta.data$LTAg_group, useNA = "always")
  # table(seurat@meta.data$LTAg_group_genotype, useNA = "always")
  
  return(seurat)
}


#' addUpperLowerBound_to_Seurat
#'
#' This function adds three columns to slot @meta.data of a Seurat object. They are 1) Filter_high_mRNA, 2) Filter_low_mRNA, and 3) Filter_ok and all three are logical (T or F) vectors. Function uses same calculation as Monocle documentation to identify cells with significantly low and high Total_mRNA counts.
#' 
#' From Monocle documentation: removed the cells with very low mRNA recovery or far more mRNA that the typical cell. Often, doublets or triplets have roughly twice the mRNA recovered as true single cells, so the latter filter is another means of excluding all but single cells from the analysis. Such filtering is handy if your protocol doesn't allow directly visualization of cell after they've been captured. 
#' @param seurat A Seurat object
#' @param pheno.use The @meta.data column name to use (default = 'Total_mRNA')
#' @keywords Seurat
#' @import Seurat
#' @return A Seurat object with 3 new columns in @meta.data
#' @export
#' @examples
#' seurat = addUpperLowerBound_to_Seurat(seurat)
#' seurat = addUpperLowerBound_to_Seurat(seurat, pheno.use = 'nUMI')

addUpperLowerBound_to_Seurat = function(seurat, pheno.use = "Total_mRNA") {
  require(Seurat)
  m = mean(log10(seurat@meta.data[,pheno.use]))
  s = sd(log10(seurat@meta.data[,pheno.use]))
  upper_bound = 10^(m + 2*s)
  lower_bound = 10^(m - 2*s)
  
  # update 'local' attributes
  attrlist = attributes(seurat)$local
  if (is.null(attrlist)) {
    attrlist = list()
  }
  attrlist[['filter.pheno.use']] = pheno.use
  attrlist[['filter.upper.bound']] = upper_bound
  attrlist[['filter.lower.bound']] = lower_bound
  attributes(seurat)$local = attrlist

  filter = data.frame("Filter_high_mRNA" = seurat@meta.data[,pheno.use] >= upper_bound,
                      "Filter_low_mRNA" = seurat@meta.data[,pheno.use] <= lower_bound)
  rownames(filter) = rownames(seurat@meta.data)
  head(filter)
  filter$Filter_ok = (filter$Filter_high_mRNA == F & filter$Filter_low_mRNA == F)
  seurat = AddMetaData(seurat, metadata = filter)
  
  return(seurat)
}


#' densityPlot_Seurat
#'
#' Outputs a Density plot for all cells. The data used will be the same metadata column name used in the addUpperLowerBound_to_Seurat() function call. Draws vertical lines to show the lower and upper bound cutoff.
#' @param seurat A Seurat object
#' @keywords Seurat
#' @import Seurat
#' @return Nothing
#' @export
#' @examples
#' densityPlot_Seurat(seurat)

densityPlot_Seurat = function(seurat) {
  require(Seurat)
  pheno.use = attributes(seurat)$local[['filter.pheno.use']]
  ub = attributes(seurat)$local[['filter.upper.bound']]
  lb = attributes(seurat)$local[['filter.lower.bound']]
  
  d = density(seurat@meta.data[,pheno.use])
  plot(d, main = pheno.use)
  abline(v = lb, col = 4)
  abline(v = ub, col = 2)
}

#' output_process_heatmaps
#'
#' Wrapper around process_heatmap() that outputs all heatmaps to a PDF.
#' @param seurat_align A Seurat alignment object
#' @param markers Marker table from Seurat 
#' @param pdffile PDF filename to save results
#' @keywords Seurat processheatmaps heatmap
#' @import gridExtra
#' @return Nothing
#' @export
#' @examples
#' output_process_heatmaps(aligned, markers, "processheatmaps.alignment.pdf")

output_process_heatmaps = function (seurat_align, markers, pdffile, ...) {
  # Create heatmaps for each Process in marker and marker.conserved tables
  #   markers
  require(gridExtra)
  while (!is.null(dev.list()))  dev.off()
  pdf(pdffile)
  processes = names(table(markers$Process))[!grepl(";", names(table(markers$Process)))]
  for (process in processes) {
    print (process)
    p = process_heatmap(seurat_align, markers = markers, process = process, ...)
    if (!is.null(p)) {
      grid.arrange(grobs = p[c(1:2)], ncol=2, top=process)
    }
  }
  dev.off()
}


#' process_heatmap
#'
#' Builds an unclustered and a clustered scaled heatmap for genes in a Process from Marker tables.
#' @param object A Seurat object
#' @param markers Marker table from Seurat 
#' @param process Process name (from gene annotations)
#' @param fontsizeRow Font size of row labels in Pheatmap
#' @keywords Seurat pheatmap
#' @import Seurat dplyr pheatmap RColorBrewer
#' @return List of 3 objects: 2 pheatmaps (unclustered and clustered) and scaled data
#' @export
#' @examples
#' process_heatmaps(seurat, markers, process = process, ...)

process_heatmap = function (object, markers, process, colors = NULL, fontsizeRow = 4, ...) {
  require(dplyr)
  require(Seurat)
  
  genes = markers %>% dplyr::filter(Process == process, !is.na(gene)) %>% select(gene) %>% unlist(use.names = F) %>% unique()
  if (length(genes) < 2) {
    return (NULL)
  }

  message("Missing genes")
  index_gene_found = genes %in% rownames(object@data)
  message(paste0(genes[!index_gene_found], collapse=", "))
  message()
  genes = genes[index_gene_found]

  message("Getting average expression for each gene")
  expr = AverageExpression(object, genes.use = genes)
  # scale(expr)  # don't want this since it scales and centers rows and columns
  message()

  message("Genes not expressed in all groups or single cells")
  rowsumzero.index = rowSums(expr) == 0
  zeroexpr.genes = rownames(expr[rowsumzero.index,])
  message(paste0(zeroexpr.genes, collapse=", "))
  expr = expr[!rowsumzero.index,]  # remove genes with rowSum expression == 0
  message()
  
  message("Genes use: ", paste0(rownames(expr), collapse=", "))

  require(pheatmap)
  require(RColorBrewer)
  if (is.null(colors)) {
    #colors = colorRampPalette(c("black","skyblue2","gold"))(100)
    colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)  # default
  }
  plot1 = pheatmap(expr, color = colors, cluster_cols = F, cluster_rows=F, scale = "row", main='Unclustered', silent=T, fontsize_row=fontsizeRow, ...)  #  built in z-score calculation in pheatmap ('scale = row')
  plot2 = pheatmap(expr, color = colors,scale = "row", main="Clustered", silent=T, fontsize_row=fontsizeRow, ...)

  # you can't extract clustering data from Pheatmap so need to do it manually
  # https://stackoverflow.com/questions/32784646/can-you-extract-the-data-matrix-from-pheatmap-in-r
  expr2 = t(scale(t(expr)))
  rowclust = hclust(dist(expr2))
  reordered = expr2[rowclust$order,]
  colclust = hclust(dist(t(expr2)))
  reordered = reordered[,colclust$order]
  reordered = tbl_df(reordered) %>% dplyr::mutate(process = process, gene = rownames(reordered)) %>% select(process, gene, everything())
  
  toReturn = list(plot1[[4]], plot2[[4]], reordered)  # https://stackoverflow.com/questions/39590849/using-a-pheatmap-in-arrangegrob  (you need slot 4 !! )
  return(toReturn)
}


# Returns list of 2 pheatmaps (unclustered and clustered)
process_heatmap_old = function (object, markers, process, colors = NULL, fontsizeRow = 4, ...) {
  genes = markers %>% filter(Process == process) %>% select(gene) %>% unlist(use.names = F) %>% unique()
  if (length(genes) < 2) {
    return (NULL)
  }
  
  expr = AverageExpression(object, genes.use = genes)
  # scale(expr)  # don't want this since it scales and centers rows and columns
  
  require(pheatmap)
  require(RColorBrewer)
  if (is.null(colors)) {
    #colors = colorRampPalette(c("black","skyblue2","gold"))(100)
    colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)  # default
  }
  plot1 = pheatmap(expr, color = colors, cluster_cols = F, cluster_rows=F, scale = "row", main='Unclustered', silent=T, fontsize_row=fontsizeRow, ...)  #  built in z-score calculation in pheatmap ('scale = row')
  plot2 = pheatmap(expr, color = colors,scale = "row", main="Clustered", silent=T, fontsize_row=fontsizeRow, ...)
  
  toReturn = list(plot1[[4]], plot2[[4]])  # https://stackoverflow.com/questions/39590849/using-a-pheatmap-in-arrangegrob  (you need slot 4 !! )
  return(toReturn)
}

