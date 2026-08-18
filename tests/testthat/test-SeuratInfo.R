test_that("SeuratInfo reports assays, reductions and idents for a v3 object", {
  expect_snapshot(SeuratInfo(pbmc_small))
})

test_that("SeuratInfo reports the same summary for a v5 object", {
  expect_snapshot(SeuratInfo(pbmc_small_v5))
})
