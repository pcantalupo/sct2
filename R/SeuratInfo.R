#' @title SeuratInfo
#' @description
#' Show information about the Seurat object such as a table of the Idents and the first two rows of metadata. In addition, shows the available Reductions and Graphs. It shows a table of information about the Assays in the object and shows which Assay is the default.
#' @param seurat A Seurat object
#' @export
#' @import Seurat
#' @importFrom SeuratObject LayerData
#' @importFrom utils head str
#' @examples
#'
#' SeuratInfo(pbmc_small)
#'
SeuratInfo = function(seurat) {
  cat("Seurat version: ", as.character(seurat@version), "\n")

  cat("\nMetadata: ")
  cat(str(seurat[[]]))

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
  slotinfo = list(); slots = c("counts", "data", "scale.data")
  assays = names(seurat@assays)
  for (assay in assays) {
    defaultassay = ifelse (DefaultAssay(seurat) == assay, "YES", "")

    assaydata = lapply(slots, function(slot) {
      if (substr(as.character(seurat@version), 1,1) == "5") {
        data = LayerData(seurat, layer = slot, assay = assay)    # version 5
      } else {
        data = GetAssayData(seurat, layer = slot, assay = assay) # previous versions
      }
      return(data)
    })
    names(assaydata) = slots
    dims = unlist(lapply(slots, \(s) { paste0(nrow(assaydata[[s]]), "x", ncol(assaydata[[s]])) }))
    hvgs = length(VariableFeatures(seurat, assay = assay))

    slotinfo[[assay]] = c(defaultassay, dims, hvgs)
  }
  df = data.frame(do.call(rbind, slotinfo))
  colnames(df) = c("default", "counts", "data", "scale.data", "HVGs")
  print(df)
}

