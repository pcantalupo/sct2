#' @title SeuratInfo
#' @description
#' Show information about the Seurat object such as a table of the Idents and the first two rows of metadata. In addition, shows the available Reductions and Graphs. It shows a table of information about the Assays in the object and shows which Assay is the default.
#' @param seurat A Seurat object
#' @export
#' @import Seurat
#' @importFrom utils head
#' @examples
#' \dontrun{
#' SeuratInfo(seurat)
#' }
SeuratInfo = function(seurat) {
  message("\nSeurat object level info")
  message("------------------------")

  message("\nSeurat version: ", seurat@version)

  message("\nMetadata: ")
  print(head(seurat[[]], n=2))

  message(paste0("\nReductions: ", paste(names(seurat@reductions), collapse = ", ")))
  message(paste0("\nGraphs: ", paste(names(seurat@graphs), collapse = ", ")))
  message("\nIdents (aka levels): ")
  tab = table(Idents(seurat))    #print(table(Idents(seurat)))   # stored in pbmc@active.ident; can also use levels(seurat)
  df = data.frame(Count = as.integer(tab))
  rownames(df) = rownames(tab)
  print(t(df))

  message("\nAssays")# (default: ", DefaultAssay(seurat), ")")
  message("------")
  slotinfo = list(); slots = c("counts", "data", "scale.data")
  assays = names(seurat@assays)
  for (assay in assays) {
    defaultassay = ifelse (DefaultAssay(seurat) == assay, "YES", "")

    assaydata = lapply(slots, \(s) { GetAssayData(seurat, slot = s, assay = assay) })
    names(assaydata) = slots
    dims = unlist(lapply(slots, \(s) { paste0(nrow(assaydata[[s]]), "x", ncol(assaydata[[s]])) }))
    hvgs = length(VariableFeatures(seurat, assay = assay))

    slotinfo[[assay]] = c(defaultassay, dims, hvgs)
  }
  df = data.frame(do.call(rbind, slotinfo))
  colnames(df) = c("default", "counts", "data", "scale.data", "HVGs")
  print(df)
}

