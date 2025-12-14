#' @title SeuratInfo
#' @description
#' Show information about the Seurat object such as a table of the Idents and the first two rows of metadata. In addition, shows the available Reductions and Graphs. It shows a table of information about the Assays in the object and shows which Assay is the default.
#' @param seurat A Seurat object
#' @param metadata Show metadata? (default FALSE)
#' @export
#' @importFrom SeuratObject LayerData GetAssayData DefaultAssay VariableFeatures
#' @importFrom utils head str
#' @examples
#'
#' SeuratInfo(pbmc_small)
#'
SeuratInfo = function(seurat, metadata = FALSE) {
  cat("Seurat version: ", as.character(seurat@version), "\n")

  # Metadata
  if (metadata) {
    cat("\nMetadata: ")
    cat(str(seurat[[]]))
  }

  # Graphs
  cat(paste0("\nGraphs: ", paste(names(seurat@graphs), collapse = ", ")))

  # Reductions
  assayused = sapply(names(seurat@reductions), function(name) {
    return(seurat[[name]]@assay.used)
  })
  cat(paste0("\nReductions: ", paste(names(assayused), " (", assayused, ")", sep = "",
                                     collapse = ", ")))

  # Images
  if (length(seurat@images) > 0) {
    cat("\nImages: ", paste(names(seurat@images), collapse = ", "))
  }

  # Idents
  cat("\nIdent label:", FindIdentLabel(seurat))
  cat("\nIdents():\n")
  tab = table(Idents(seurat))
  df = data.frame(Count = as.integer(tab))
  rownames(df) = rownames(tab)
  print(t(df))

  # Assays
  cat("\nAssays:\n")
  assays = names(seurat@assays)
  default_assay = DefaultAssay(seurat)

  # get the dimension string for a slot in an assay object (i.e. "6182x198075")
  get_layer_dim = function(assay_obj, slot) {
    slotnames = slotNames(assay_obj)
    layer = NULL
    if ("layers" %in% slotnames) {   # v5 layers slot
      layer = assay_obj@layers[[slot]]
    }
    else if (slot %in% slotnames) {  # v4 direct slot
      layer = slot(assay_obj, slot)
    }

    if (is.null(layer)) {
      dim_string = "0x0"
    } else {
      dims = dim(layer)
      dim_string = paste0(dims[1], "x", dims[2])
    }

    return(dim_string)
  }

  slotinfo = list()
  slots = c("counts", "data", "scale.data")
  assay = "RNA"
  for (assay in assays) {
    assay_obj = seurat@assays[[assay]]  # get the Assay object

    defaultassay = ifelse (default_assay == assay, "YES", "")

    dims = sapply(slots, function(slot) get_layer_dim(assay_obj, slot))

    hvgs = length(VariableFeatures(seurat, assay = assay))

    slotinfo[[assay]] = c(defaultassay, dims, hvgs)
  }

  df = data.frame(do.call(rbind, slotinfo))
  colnames(df) = c("default", slots, "HVGs")
  print(df)
}

