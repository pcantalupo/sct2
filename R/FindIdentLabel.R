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
#'
#' # Get the identity label for pbmc_small
#' FindIdentLabel(pbmc_small)
#'
FindIdentLabel <- function(seurat) {
  ident.label <- as.character(Idents(seurat))
  chr_cols <- colnames(seurat@meta.data)[
    sapply(seurat@meta.data, function(x) is.character(x) || is.factor(x))
  ]
  labels <- sapply(chr_cols,
                   function(x) all(ident.label == as.character(seurat@meta.data[, x]))) %>% .[.] %>% .[!is.na(.)]
  label <- names(labels[labels])
  label = label[!(label %in% c("seurat_clusters", "ident"))][1]
  return(label)
}

