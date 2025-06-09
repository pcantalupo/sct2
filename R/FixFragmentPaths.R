#' @title FixFragmentPaths
#' @description
#' Fix the paths to the atac_fragments.tsv.gz file in a Seurat object 
#' @param seurat A Seurat object
#' @param cellranger_path Path to the cellranger folder that contains subfolders matching the sample names in the Seurat metadata `orig.ident`
#' @param assay The assay name for the ChromatinAssay in the Seurat object
#' @param sample_col The column name in Seurat metadata that contains the sample identifiers
#' @return A Seurat object
#' @export
#' @import Seurat
#' @import Signac
#' @examples
#' \dontrun{
#' FixFragmentPaths(seurat, cellranger_path = "../cellranger_count")
#' }
#'
FixFragmentPaths <- function(seurat, cellranger_path, assay = "ATAC", sample_col = "orig.ident") {
  # Parameter checks
  if (!dir.exists(cellranger_path)) {
    stop(paste("Cellranger path", cellranger_path, "does not exist"))
  }
  if (!sample_col %in% colnames(seurat@meta.data)) {
    stop(paste("Sample column", sample_col, "not found in Seurat metadata"))
  }
  if (!assay %in% Assays(seurat)) {
    stop(paste("Assay", assay, "not found in Seurat object"))
  }

  # Retrieve the ChromatinAssay
  chromassay <- seurat[[assay]]
  frags <- Fragments(chromassay)
  Fragments(chromassay) <- NULL  # Temporarily remove fragments
  
  # Process each fragment
  for (i in seq_along(frags)) {
    frag <- frags[[i]]
    cells_in_frag <- GetFragmentData(frag, slot = "cells")
    cell_names <- names(cells_in_frag)

    # Some fragments might be empty if Seurat object was subsetted on a particular sample
    if (length(cell_names) == 0) {
      warning("There are no cells in fragment ", i, "...skipping this fragment")
      next
    } 

    # Get sample from metadata for the fragment's cells
    samples <- seurat@meta.data[cell_names, sample_col]
    unique_samples <- unique(samples)
    
    # Ensure all cells in the fragment belong to one sample
    if (length(unique_samples) != 1) {
      stop(paste("Cells in fragment", i, "belong to multiple samples:", paste(unique_samples, collapse = ", ")))
    }
    sample <- unique_samples[1]
    message("\nSample: ", sample)
    
    # Define current and new fragment paths
    current_frag_path <- frag@path
    new_frag_path <- normalizePath(file.path(cellranger_path, sample, "outs/atac_fragments.tsv.gz"), mustWork = FALSE)
    
    # Update path if it doesn’t exist or is incorrect
    if (!file.exists(current_frag_path) || current_frag_path != new_frag_path) {
      if (!file.exists(new_frag_path)) {
        stop(paste("New fragment path does not exist:", new_frag_path))
      }
      message("Updating fragment path from: ", current_frag_path)
      message("to: ", new_frag_path)
      frag <- UpdatePath(frag, new.path = new_frag_path)
    } else {
      message("Fragment path is already correct: ", current_frag_path)
    }
    
    frags[[i]] <- frag
  }
  
  # Update the ChromatinAssay and Seurat object
  Fragments(chromassay) <- frags
  seurat[[assay]] <- chromassay
  return(seurat)
}
