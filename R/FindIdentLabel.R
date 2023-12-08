#' @title FindIdentLabel
#' @description
#' Find identical label between ident and metadata
#' @references https://github.com/nyuhuyang/SeuratExtra/blob/master/R/Seurat4_functions.R
#' @param seurat A Seurat object
#' @return The colname in metadata
#' @export
#' @importFrom Seurat Idents 
#' @importFrom magrittr %>%
#' @examples
#' \dontrun{
#' FindIdentLabel(seurat)
#' }
FindIdentLabel <- function(seurat) {
  ident.label <- as.character(Idents(seurat))
  labels <- sapply(colnames(seurat@meta.data),
                   function(x) all(ident.label == seurat@meta.data[,x])) %>% .[.] %>% .[!is.na(.)]
  label <- names(labels[labels])
  label = label[!(label %in% c("seurat_clusters","ident"))][1]
  return(label)
}

