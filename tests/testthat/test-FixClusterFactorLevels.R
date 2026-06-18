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
