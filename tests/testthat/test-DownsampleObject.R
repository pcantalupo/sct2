
test_that("DownsampleObject returns the exact same object when downsample = 1", {
  expect_equal(pbmc_small, DownsampleObject(pbmc_small, 1))
  expect_equal(pbmc_small_sce, DownsampleObject(pbmc_small_sce, 1))
})

test_that("a value of 1 or less is a fraction of the cells", {
  expect_equal(ncol(DownsampleObject(pbmc_small, 0.1)), 8)
  expect_equal(ncol(DownsampleObject(pbmc_small_sce, 0.25)), 20)
  expect_equal(ncol(DownsampleObject(multiome_small, 0.5)), 60)
})

test_that("a fractional number of cells is rounded down", {
  # 80 * 0.33 = 26.4
  expect_equal(ncol(DownsampleObject(pbmc_small, 0.33)), 26)
})

test_that("a value above 1 is an absolute number of cells", {
  expect_equal(ncol(DownsampleObject(pbmc_small, 20)), 20)
  expect_equal(ncol(DownsampleObject(pbmc_small_sce, 20)), 20)
})

test_that("2 is the smallest absolute number of cells", {
  expect_equal(ncol(DownsampleObject(pbmc_small, 2)), 2)
})

test_that("an absolute number at or above the object size returns the same object", {
  expect_equal(pbmc_small, DownsampleObject(pbmc_small, ncol(pbmc_small)))
  expect_equal(pbmc_small, DownsampleObject(pbmc_small, 500))
  expect_equal(pbmc_small_sce, DownsampleObject(pbmc_small_sce, 500))
})

test_that("the kept cells are cells of the input object", {
  result = DownsampleObject(pbmc_small, 20)
  expect_true(all(colnames(result) %in% colnames(pbmc_small)))
  expect_false(anyDuplicated(colnames(result)) > 0)
})

test_that("the same seed keeps the same cells and a different seed does not", {
  first  = DownsampleObject(pbmc_small, 20, seed = 1946)
  second = DownsampleObject(pbmc_small, 20, seed = 1946)
  other  = DownsampleObject(pbmc_small, 20, seed = 42)
  expect_equal(colnames(first), colnames(second))
  expect_false(identical(colnames(first), colnames(other)))
})

test_that("the RNG is not consumed when there is no downsampling", {
  set.seed(99)
  expected = runif(1)

  set.seed(99)
  DownsampleObject(pbmc_small)
  expect_equal(runif(1), expected)
})

test_that("a downsampled Seurat object keeps its reductions and metadata columns", {
  result = DownsampleObject(pbmc_small, 20)
  expect_s4_class(result, "Seurat")
  expect_equal(Seurat::Reductions(result), Seurat::Reductions(pbmc_small))
  expect_equal(colnames(result[[]]), colnames(pbmc_small[[]]))
  expect_equal(nrow(result), nrow(pbmc_small))
})

test_that("a downsampled SingleCellExperiment keeps its colData columns", {
  result = DownsampleObject(pbmc_small_sce, 20)
  expect_s4_class(result, "SingleCellExperiment")
  expect_equal(colnames(SummarizedExperiment::colData(result)),
               colnames(SummarizedExperiment::colData(pbmc_small_sce)))
  expect_equal(nrow(result), nrow(pbmc_small_sce))
})

test_that("a downsampled multiome object keeps both assays", {
  result = DownsampleObject(multiome_small, 20)
  expect_equal(Seurat::Assays(result), Seurat::Assays(multiome_small))
})

test_that("the metadata of the kept cells is carried over unchanged", {
  result = DownsampleObject(pbmc_small, 20)
  expect_equal(result[[]], pbmc_small[[]][colnames(result), ])
})
