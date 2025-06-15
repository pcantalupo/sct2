#' SaveMetadata
#' @description
#' Save the Seurat metadata to a TSV file.
#'
#' @param seurat A Seurat object
#' @param filename A TSV filename
#' @param colname_for_rows The name used for the column name of the cell rownames
#'
#' @export
#' @importFrom tibble rownames_to_column
#' @importFrom magrittr %>%
#' @importFrom readr write_tsv
#'
#' @examples
#' \dontrun{
#' SaveMetadata(pbmc_small, "pbmc_small_metadata.tsv")
#' }
#'
SaveMetadata = function(seurat, filename, colname_for_rows = "cellid") {
  if(missing(seurat) | missing(filename)) {
    stop("seurat and filename parameters are required.", call. = FALSE)
  }

  toWrite = seurat[[]] %>% tibble::rownames_to_column(var = colname_for_rows)
  write_tsv(toWrite, filename)
}
