test_that("FindIdentLabel works", {
  expect_equal(FindIdentLabel(pbmc_small), "RNA_snn_res.1")
})

test_that("FindIdentLabel works when metadata contains POSIXct column", {
  seurat = pbmc_small
  seurat@meta.data$datetime_col = as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  expect_equal(FindIdentLabel(seurat), "RNA_snn_res.1")
})
