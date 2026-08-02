#' SaveMetadata
#' @description
#' Save the Seurat metadata to a TSV, RDS or QS2 file. The output format is
#' chosen from the filename extension: `.rds` and `.qs2` save the metadata
#' data.frame directly (preserving factor levels, column classes and the cell
#' barcodes as rownames), any other extension writes a TSV.
#'
#' @param seurat A Seurat object
#' @param filename Output filename (.tsv, .rds or .qs2)
#' @param colname_for_rows The name used for the column name of the cell
#'   rownames. TSV output only; RDS and QS2 keep the cell barcodes as rownames.
#'
#' @return Invisibly returns filename
#'
#' @export
#' @importFrom tibble rownames_to_column
#' @importFrom magrittr %>%
#' @importFrom readr write_tsv
#' @importFrom tools file_ext
#'
#' @examples
#' \dontrun{
#' SaveMetadata(pbmc_small, "pbmc_small_metadata.tsv")
#' SaveMetadata(pbmc_small, "pbmc_small_metadata.rds")
#' SaveMetadata(pbmc_small, "pbmc_small_metadata.qs2")
#' }
#'
SaveMetadata = function(seurat, filename, colname_for_rows = "cellid") {
  if(missing(seurat) | missing(filename)) {
    stop("seurat and filename parameters are required.", call. = FALSE)
  }

  metadata = seurat[[]]
  ext = tolower(tools::file_ext(filename))
  if (ext == "rds") {
    saveRDS(metadata, filename)
  } else if (ext == "qs2") {
    qs2::qs_save(metadata, filename)
  } else {
    toWrite = metadata %>% tibble::rownames_to_column(var = colname_for_rows)
    write_tsv(toWrite, filename)
  }
  invisible(filename)
}
