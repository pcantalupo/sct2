#' @title SeuratInfo
#' @description
#' Show information about the Seurat object such as a table of the Idents and the first two rows of metadata. In addition, shows the available Reductions and Graphs. It shows a table of information about the Assays in the object and shows which Assay is the default.
#' @param seurat A Seurat object
#' @param metadata Show metadata? (default FALSE)
#' @export
#' @import Seurat
#' @importFrom SeuratObject LayerData
#' @importFrom utils head str
#' @examples
#'
#' SeuratInfo(pbmc_small)
#'
SeuratInfo = function(seurat, metadata = FALSE) {
  cat("Seurat version: ", as.character(seurat@version), "\n")

  if (metadata) {
    cat("\nMetadata: ")
    cat(str(seurat[[]]))
  }

  cat(paste0("\nGraphs: ", paste(names(seurat@graphs), collapse = ", ")))

  assayused = sapply(names(seurat@reductions), function(name) {
    return(seurat[[name]]@assay.used)
  })
  cat(paste0("\nReductions: ", paste(names(assayused), " (", assayused, ")", sep = "",
                                     collapse = ", ")))

  if (length(seurat@images) > 0) {
    cat("\nImages: ", paste(names(seurat@images), collapse = ", "))
  }

  cat("\nIdent label:", FindIdentLabel(seurat))

  cat("\nIdents():\n")
  tab = table(Idents(seurat))    #print(table(Idents(seurat)))   # stored in pbmc@active.ident; can also use levels(seurat)
  df = data.frame(Count = as.integer(tab))
  rownames(df) = rownames(tab)
  print(t(df))


  cat("\nAssays:\n")# (default: ", DefaultAssay(seurat), ")")
  is_v5 = substr(as.character(seurat@version), 1, 1) == "5"
  default_assay = DefaultAssay(seurat)

  slotinfo = list(); slots = c("counts", "data", "scale.data")
  chromatin_info = list()

  assays = names(seurat@assays)
  for (assay in assays) {
    defaultassay = ifelse (default_assay == assay, "YES", "")

    # Get dimensions efficiently without returning full matrices
    dims = sapply(slots, function(slot) {
      tryCatch({
        if (is_v5) {
          layer_data = LayerData(seurat, layer = slot, assay = assay)
        } else {
          layer_data = GetAssayData(seurat, layer = slot, assay = assay)
        }
        paste0(nrow(layer_data), "x", ncol(layer_data))
      }, error = function(e) {
        "0x0"  # Handle missing slots gracefully
      })
    })

    hvgs = length(VariableFeatures(seurat, assay = assay))

    slotinfo[[assay]] = c(defaultassay, dims, hvgs)
  }

  df = data.frame(do.call(rbind, slotinfo))
  colnames(df) = c("default", "counts", "data", "scale.data", "HVGs")
  print(df)
}

