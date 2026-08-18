#' @title SeuratInfo
#' @description
#' Show information about the Seurat object such as a table of the Idents and the first two rows of metadata. In addition, shows the available Reductions and Graphs. It shows a table of information about the Assays in the object and shows which Assay is the default.
#'
#' The report is printed and the gathered information is returned invisibly as a
#' `seurat_info` object, so the printing behaviour is unchanged at every call site.
#' @param seurat A Seurat object
#' @param metadata Show metadata? (default FALSE)
#' @return Invisibly, an object of class `seurat_info`: a list with elements
#'   `version` (character), `metadata` (the metadata data.frame when
#'   `metadata = TRUE`, otherwise `NULL`), `graphs` (character), `reductions`
#'   (character, named by reduction name, values are the assay used), `images`
#'   (character), `ident_label` (character), `idents` (a table of cell counts per
#'   ident) and `assays` (a data.frame with one row per assay). Formatting lives
#'   in [print.seurat_info()].
#' @export
#' @importFrom SeuratObject LayerData GetAssayData DefaultAssay VariableFeatures
#' @importFrom utils head str
#' @examples
#'
#' SeuratInfo(pbmc_small)
#'
SeuratInfo = function(seurat, metadata = FALSE) {
  version = as.character(seurat@version)

  # Metadata
  if (metadata) {
    metadata_df = seurat[[]]
  } else {
    metadata_df = NULL
  }

  # Graphs
  graphs = names(seurat@graphs)

  # Reductions
  assayused = vapply(names(seurat@reductions), function(name) {
    return(seurat[[name]]@assay.used)
  }, character(1))

  # Images
  images = names(seurat@images)

  # Idents
  ident_label = FindIdentLabel(seurat)
  idents = table(Idents(seurat))

  # Assays
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

  x = structure(list(version = version,
                     metadata = metadata_df,
                     graphs = graphs,
                     reductions = assayused,
                     images = images,
                     ident_label = ident_label,
                     idents = idents,
                     assays = df),
                class = "seurat_info")

  print(x)
  invisible(x)
}
