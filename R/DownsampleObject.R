#' @title DownsampleObject
#' @description
#' Draw a uniform random subset of cells from a single cell object. Works on any
#' object with cells as columns, including Seurat and SingleCellExperiment, since
#' the subset is done with `object[, cells]`.
#'
#' The `downsample` argument is overloaded: a value of 1 or less is a fraction of
#' the cells to keep, and a value above 1 is an absolute cell count (so the
#' smallest count is 2). A count larger than the object is clamped to the object
#' size. A fraction small enough to resolve to zero cells is an error.
#'
#' When the resolved number of cells to keep is the whole object -- the default
#' `downsample = 1`, or a count at or above the object size -- the object is
#' returned unchanged and the RNG is not used.
#' @param object A Seurat or SingleCellExperiment object
#' @param downsample Fraction of cells to keep (<= 1) or absolute number of cells
#'   to keep (> 1) [default: 1, no downsampling]
#' @param seed RNG seed for the sample
#' @return An object of the same class as `object`, subset to the sampled cells
#' @export
#' @examples
#'
#' DownsampleObject(pbmc_small, 0.1)
#' DownsampleObject(pbmc_small_sce, 20)
#'
DownsampleObject = function(object, downsample = 1, seed = 1946) {
  n_cells = ncol(object)
  if (downsample <= 1) {
    n_keep = floor(n_cells * downsample)
  } else {
    n_keep = min(downsample, n_cells)
  }

  if (n_keep == n_cells) {
    return(object)
  }
  if (n_keep < 1) {
    stop("Resolved cells-to-keep is ", n_keep, "; nothing to sample")
  }

  set.seed(seed)
  cells = sample(colnames(object), size = n_keep, replace = FALSE)
  object_ds = object[, cells]

  return(object_ds)
}
