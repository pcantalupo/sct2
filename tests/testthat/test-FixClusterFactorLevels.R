test_that("FixClusterFactorLevels relevels factor columns in numerical order", {
  # multiome_small ships LexSCT_snn_res.0.8 with its 18 levels in lexicographic
  # order, so a function that returns the object untouched fails this test
  seurat = multiome_small
  expect_equal(levels(seurat$LexSCT_snn_res.0.8),
               c("0", "1", "10", "11", "12", "13", "14", "15", "16", "17",
                 "2", "3", "4", "5", "6", "7", "8", "9"))

  result = suppressMessages(FixClusterFactorLevels(seurat))

  expect_equal(levels(result$LexSCT_snn_res.0.8),
               c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
                 "10", "11", "12", "13", "14", "15", "16", "17"))

  # the remaining resolution columns arrive sorted and must stay that way
  res_cols = grep("snn_res", colnames(result@meta.data), value = TRUE)
  expect_length(res_cols, 14)
  for (col in res_cols) {
    lvls = levels(result@meta.data[[col]])
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
  # diverge. multiome_small carries LexSCT_snn_res.0.8, the same 18 cluster
  # assignments as SCT_snn_res.0.8 but with levels in lexicographic order.
  seurat = multiome_small
  expected = seurat$SCT_snn_res.0.8

  # confirm the fixture actually reproduces the problem before testing the fix
  expect_false(identical(levels(seurat$LexSCT_snn_res.0.8), levels(expected)))

  result = suppressMessages(FixClusterFactorLevels(seurat))

  expect_equal(levels(result$LexSCT_snn_res.0.8), levels(expected))
  # releveling must reorder levels only, never reassign a cell
  expect_equal(as.character(result$LexSCT_snn_res.0.8), as.character(expected))
})
