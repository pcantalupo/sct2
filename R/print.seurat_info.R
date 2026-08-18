#' @title Print a seurat_info object
#' @description
#' Formats and prints the report gathered by [SeuratInfo()]. All of the
#' formatting for that report lives here; `SeuratInfo()` only gathers the data.
#' @param x A `seurat_info` object, as returned by [SeuratInfo()]
#' @param ... Ignored
#' @return `x`, invisibly
#' @export
#' @importFrom utils str
print.seurat_info = function(x, ...) {
  cat("Seurat version: ", x$version, "\n")

  # Metadata
  if (!is.null(x$metadata)) {
    cat("\nMetadata: ")
    str(x$metadata)
  }

  # Graphs
  cat(paste0("\nGraphs: ", paste(x$graphs, collapse = ", ")))

  # Reductions
  cat(paste0("\nReductions: ", paste(names(x$reductions), " (", x$reductions, ")", sep = "",
                                     collapse = ", ")))

  # Images
  if (length(x$images) > 0) {
    cat("\nImages: ", paste(x$images, collapse = ", "))
  }

  # Idents
  cat("\nIdent label:", x$ident_label)
  cat("\nIdents():\n")
  df = data.frame(Count = as.integer(x$idents))
  rownames(df) = rownames(x$idents)
  print(t(df))

  # Assays
  cat("\nAssays:\n")
  print(x$assays)

  invisible(x)
}
