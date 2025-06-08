#' PBMC Seurat object
#'
#' This is the pbmc_small Seurat object from the Seurat package
#'
#' @format
#' Seurat object with 230 rows/genes and 80 columns/cells
#' \describe{
#'   \item{nCount_RNA}{Number of UMI counts in the cell}
#'   \item{nFeature_RNA}{Number of genes detected in the cell}
#'   \item{RNA_snn_res.0.8}{Cluster identity using resolution 0.8}
#'   \item{RNA_snn_res.1}{Cluster identity using resolution 1}
#'   \item{letter.idents}{Sample id?}
#'   \item{groups}{Group id}
#' }
#' @source https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k
#' @references Seurat package https://satijalab.org/seurat/
"pbmc_small"

#' PBMC SCE object
#'
#' This is the pbmc_small Seurat object from the Seurat package converted
#' to a SingleCellExperiment object
#'
#' @format
#' Seurat object with 230 rows/genes and 80 columns/cells
#' \describe{
#'   \item{nCount_RNA}{Number of UMI counts in the cell}
#'   \item{nFeature_RNA}{Number of genes detected in the cell}
#'   \item{RNA_snn_res.0.8}{Cluster identity using resolution 0.8}
#'   \item{RNA_snn_res.1}{Cluster identity using resolution 1}
#'   \item{letter.idents}{Sample id?}
#'   \item{groups}{Group id}
#' }
#' @source https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k
#' @references Seurat package https://satijalab.org/seurat/
"pbmc_small_sce"


#' PBMC Seurat object (version 5)
#'
#' Converted pbmc_small into a version 5 Seurat object.
#'
#' @format
#' Seurat object with 230 rows/genes and 80 columns/cells.
#' \describe{
#'   \item{nCount_RNA}{Number of UMI counts in the cell}
#'   \item{nFeature_RNA}{Number of genes detected in the cell}
#'   \item{RNA_snn_res.0.8}{Cluster identity using resolution 0.8}
#'   \item{RNA_snn_res.1}{Cluster identity using resolution 1}
#'   \item{letter.idents}{Sample id?}
#'   \item{groups}{Group id}
#' }
#' @source https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k
#' @references Seurat package https://satijalab.org/seurat/
"pbmc_small_v5"


#' Multiome Seurat object
#'
#' Multiome Seurat object with 6 samples (20 cells each)
#'
#' @format
#' Contains two assays, RNA and ATAC. ATAC assay is a ChromatinAssay created with Signac
"multiome_small"

