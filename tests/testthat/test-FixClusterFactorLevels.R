test_that("FixClusterFactorLevels relevels factor columns in numerical order", {
  seurat <- pbmc_small
  result <- FixClusterFactorLevels(seurat)
  res_cols <- grep("snn_res", colnames(result@meta.data), value = TRUE)
  for (col in res_cols) {
    lvls <- levels(result@meta.data[[col]])
    expect_equal(lvls, as.character(sort(as.numeric(lvls))))
  }
})

test_that("FixClusterFactorLevels ignores numeric score columns", {
  seurat <- pbmc_small
  seurat@meta.data$RNA_snn_res.0.2 <- factor(seurat@meta.data$RNA_snn_res.1)
  seurat@meta.data$RNA_snn_res.0.2.score <- runif(ncol(seurat))
  expect_no_error(FixClusterFactorLevels(seurat))
})

test_that("FixClusterFactorLevels ignores character score columns", {
  seurat <- pbmc_small
  seurat@meta.data$RNA_snn_res.0.5.score <- as.character(runif(ncol(seurat)))
  expect_no_error(FixClusterFactorLevels(seurat))
})

test_that("FixClusterFactorLevels sorts levels that are out of numeric order", {
  # pbmc_small tops out at 3 clusters, where string and numeric order cannot
  # diverge. multiome_small has 18 clusters at this resolution, and levels go
  # out of order whenever the cluster ids pass through a character vector.
  seurat = multiome_small
  orig = seurat$SCT_snn_res.0.8
  seurat$SCT_snn_res.0.8 = factor(as.character(orig))

  # confirm the fixture actually reproduces the problem before testing the fix
  expect_false(identical(levels(seurat$SCT_snn_res.0.8), levels(orig)))

  result = suppressMessages(FixClusterFactorLevels(seurat))

  expect_equal(levels(result$SCT_snn_res.0.8), levels(orig))
  # releveling must reorder levels only, never reassign a cell
  expect_equal(as.character(result$SCT_snn_res.0.8), as.character(orig))
})
