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


  # Assays
  cat("\nAssays:\n")
  assays = names(seurat@assays)
  default_assay = DefaultAssay(seurat)
  
  # get the dimension string for a slot in an assay object (i.e. "6182x198075")
  get_layer_dim = function(asy, slot_name) {
    tryCatch({
      sn = slotNames(asy)
      lyr = NULL
      
      # v5: layers slot
      if ("layers" %in% sn && !is.null(asy@layers[[slot_name]])) {
        lyr = asy@layers[[slot_name]]
      }
      # v3/v4: direct slot
      else if (slot_name %in% sn) {
        lyr = slot(asy, slot_name)
      }
      
      if (is.null(lyr)) return("0x0")
      
      # dim() if possible
      d = tryCatch(dim(lyr), error = function(e) NULL)
      if (!is.null(d) && length(d) == 2) {
        return(paste0(d[1], "x", d[2]))
      }
      
      # SCE fallback
      if (inherits(lyr, "SingleCellExperiment")) {
        a = SummarizedExperiment::assayNames(lyr)
        use = if ("counts" %in% a) "counts" else a[1]
        d = dim(SummarizedExperiment::assay(lyr, use))
        return(paste0(d[1], "x", d[2]))
      }
      
      "unknown"
    }, error = function(e) {
      paste0("ERROR:", gsub("\n"," ", conditionMessage(e)))
    })
  }
  
  slotinfo = list(); slots = c("counts", "data", "scale.data")
  assay = "RNA"
  for (assay in assays) {
    asy = seurat@assays[[assay]]  # get the Assay object
    sn = slotNames(asy)           # get the slot names for this object
    
    defaultassay = ifelse (default_assay == assay, "YES", "")
    
    dims = sapply(slots, function(slot) get_layer_dim(asy, slot))
    
    hvgs = length(VariableFeatures(seurat, assay = assay))
    
    slotinfo[[assay]] = c(defaultassay, dims, hvgs)
  }
  
  df = data.frame(do.call(rbind, slotinfo))
  colnames(df) = c("default", slots, "HVGs")
  print(df)
}

