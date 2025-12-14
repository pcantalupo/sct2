summarize_seurat <- function(seurat, metadata = FALSE) {

  cat("Seurat version:", as.character(seurat@version), "\n")
  
  if (metadata) {
    cat("\nMetadata:\n")
    print(str(seurat[[]]))
  }
  
  # Graphs
  cat("\nGraphs:", paste(names(seurat@graphs), collapse = ", "))
  
  # Reductions (safe)
  reduction_assays <- sapply(names(seurat@reductions), function(rname) {
    red <- seurat@reductions[[rname]]
    if (!is.null(red)) {
      if ("assay.used" %in% slotNames(red)) {
        return(red@assay.used %||% NA_character_)
      }
    }
    NA_character_
  })

  cat("\nReductions:",
      paste(names(reduction_assays), " (", reduction_assays, ")",
            sep = "", collapse = ", "))
  
  # Images
  if (length(seurat@images) > 0)
    cat("\nImages:", paste(names(seurat@images), collapse = ", "))
  
  # Idents
  cat("\nIdent label:", FindIdentLabel(seurat))
  cat("\nIdents():\n")
  tab <- table(Idents(seurat))
  df <- data.frame(Count = as.integer(tab))
  rownames(df) <- names(tab)
  print(t(df))
  
  # Assays
  cat("\nAssays:\n")
  assays = names(seurat@assays)
  default_assay = DefaultAssay(seurat)
  
  # get the dimension string for a slot in an assay object (i.e. "6182x198075")
  get_layer_dim <- function(asy, slot_name) {
    tryCatch({
      sn <- slotNames(asy)
      lyr <- NULL
      
      # v5: layers slot
      if ("layers" %in% sn && !is.null(asy@layers[[slot_name]])) {
        lyr <- asy@layers[[slot_name]]
      }
      # v3/v4: direct slot
      else if (slot_name %in% sn) {
        lyr <- slot(asy, slot_name)
      }
      
      if (is.null(lyr)) return("0x0")
      
      # dim() if possible
      d <- tryCatch(dim(lyr), error = function(e) NULL)
      if (!is.null(d) && length(d) == 2) {
        return(paste0(d[1], "x", d[2]))
      }
      
      # SCE fallback
      if (inherits(lyr, "SingleCellExperiment")) {
        a <- SummarizedExperiment::assayNames(lyr)
        use <- if ("counts" %in% a) "counts" else a[1]
        d <- dim(SummarizedExperiment::assay(lyr, use))
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
    asy <- seurat@assays[[assay]]  # get the Assay object
    sn <- slotNames(asy)           # get the slot names for this object
    
    defaultassay = ifelse (default_assay == assay, "YES", "")
    
    dims <- sapply(slots, function(slot) get_layer_dim(asy, slot))
    
    hvgs = length(VariableFeatures(seurat, assay = assay))
    
    slotinfo[[assay]] <- c(defaultassay, dims, hvgs)
  }
  
  df <- data.frame(do.call(rbind, slotinfo))
  colnames(df) <- c("default", slots, "HVGs")
  print(df)
}
