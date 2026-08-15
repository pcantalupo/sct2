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
#' 10x Multiome object with 6 samples (20 cells each), subsampled to 100 RNA
#' features and 100 ATAC peaks. `Idents()` is `orig.ident`.
#'
#' @format
#' Seurat object with 200 features across 120 cells in two assays.
#'
#' Assays:
#' \describe{
#'   \item{RNA}{`Assay5`, the default assay; 100 features; layers `counts`, `data`, `scale.data`}
#'   \item{ATAC}{Signac `ChromatinAssay`; 100 peaks named `chrom-start-end`; gene
#'     annotation present; 6 `Fragment` objects whose paths are relative
#'     (`cellranger_count/<sample>/outs/atac_fragments.tsv.gz`) and do not exist on disk}
#' }
#' The SCT assay has been dropped, but the `nCount_SCT`/`nFeature_SCT` columns,
#' the SCT-derived cluster columns and both reductions still reference it.
#'
#' Dimensional reductions (2 dims each, both with assay `SCT`):
#' \describe{
#'   \item{SCT_pca_umap}{UMAP on `SCT_pca`; key `SCTpcaUMAP_`}
#'   \item{umap}{UMAP on `harmony`; key `umap_`}
#' }
#'
#' Sample and design metadata:
#' \describe{
#'   \item{orig.ident}{Sample id; factor with levels `2`, `3`, `4`, `5a`, `5b`, `6`; 20 cells each}
#'   \item{RetinalArea}{`VT+DT`, `EquatorN`, `EquatorT`, `HAA`}
#'   \item{EmbryoStage}{`HH26/27` (100 cells), `HH29` (20 cells)}
#' }
#'
#' QC metadata:
#' \describe{
#'   \item{nCount_RNA, nFeature_RNA, nCount_SCT, nFeature_SCT}{RNA counts and genes detected}
#'   \item{percent.mt, log10GenesPerUMI}{RNA QC metrics}
#'   \item{nCount_ATAC, nFeature_ATAC, atac_barcode, atac_fragments, atac_TSS_fragments,
#'     atac_peak_region_fragments, atac_peak_region_cutsites}{Cell Ranger ARC ATAC metrics}
#'   \item{nucleosome_signal, nucleosome_percentile, TSS.enrichment, TSS.percentile,
#'     pct_reads_in_peaks, blacklist_ratio}{Signac ATAC QC metrics}
#'   \item{scDblFinder.class}{RNA doublet call; all 120 cells are `singlet`}
#'   \item{AMULET_class}{ATAC doublet call; 119 `singlet` and 1 `NA`, the only
#'     missing value in the object}
#' }
#'
#' Clustering metadata. Every column is a factor, and all of them except
#' `LexSCT_snn_res.0.8` have their levels in numerical order
#' (`0`, `1`, ..., `n-1`) rather than lexicographic order:
#' \describe{
#'   \item{SCT_snn_res.0.1 - SCT_snn_res.0.8}{Clusters on the SCT PCA graph;
#'     7, 9, 12, 13, 14, 15, 17 and 18 clusters at resolutions 0.1 - 0.8}
#'   \item{HarmSCT_snn_res.0.1 - HarmSCT_snn_res.0.5}{Clusters on the harmony graph;
#'     5, 8, 9, 12 and 13 clusters at resolutions 0.1 - 0.5}
#'   \item{seurat_clusters}{13 clusters; the same assignments as `HarmSCT_snn_res.0.5`}
#'   \item{LexSCT_snn_res.0.8}{The same assignments as `SCT_snn_res.0.8`, but with
#'     the 18 levels in lexicographic order (`0`, `1`, `10`, `11`, ..., `17`, `2`, ..., `9`).
#'     Added so cluster level ordering can be tested}
#' }
#' No `Graphs()` are stored. `@commands` retains the 11 Seurat calls used to build the object.
"multiome_small"

