library(Seurat)
test_that("FindIdentLabel works", {
  expect_equal(FindIdentLabel(pbmc_small), "RNA_snn_res.1")
})
