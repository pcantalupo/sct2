#' SignacInfo
#' @description
#' Show information about the Signac ChromatinAssays within a Seurat object
#' @param seurat A Seurat object
#' @export
#' @import Seurat
#' @importFrom GenomicRanges width
#' @importFrom utils head str
#'
#' @examples
#' \dontrun{
#' SignacInfo(pbmc_small)
#' }
SignacInfo = function(seurat) {
  chromatin_info = list()
  assays = Assays(seurat)

  for (assay in assays) {
    assay_obj = seurat@assays[[assay]]
    if (inherits(assay_obj, "ChromatinAssay")) {

      # Peak width statistics
      peak_widths = width(granges(assay_obj))
      width_summary = summary(peak_widths)

      # Fragment information
      fragments = Fragments(assay_obj)
      frag_info = lapply(fragments, function(frag) {
        ncells = length(frag@cells)
        path = frag@path
        list(ncells = ncells, path = path)
      })

      chromatin_info[[assay]] = list(
        peaks = nrow(assay_obj),
        cells = ncol(assay_obj),
        width_min = width_summary[1],
        width_median = width_summary[3],
        width_max = width_summary[6],
        fragments = length(fragments),
        frag_info = frag_info)
    }
  }

  # Check if any ChromatinAssays found
  if (length(chromatin_info) == 0) {
    cat("No ChromatinAssays found in this Seurat object.\n")
    return(invisible())
  }

  # Print ChromatinAssay info
  cat("\nChromatin Assays:\n\n")
  for (assay_name in names(chromatin_info)) {
    info = chromatin_info[[assay_name]]

    info_to_print = paste0(assay_name, ": ", info$peaks, " peaks, ", info$cells, " cells,",
                           " width range [", info$width_min, "-", info$width_max, "],",
                           " median ", info$width_median, ", fragments: ", info$fragments)
    cat(info_to_print, "\n", sep = "")

    for(i in 1:length(info$frag_info)) {
      cat(paste0("  ", info$frag_info[[i]]$path, " (", info$frag_info[[i]]$ncells,
                 " cells)", "\n"))
    }
    cat("\n")
  }

  invisible(NULL)
}

