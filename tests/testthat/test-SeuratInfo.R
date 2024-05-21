test_that("SeuratInfo works", {
  #val = SeuratInfo(pbmc_small)
  #dput(val)
  # checking that the last print() in SeuratInfo returns the correct value
  expect_equal(SeuratInfo(pbmc_small),
               structure(list(default = "YES", counts = "230x80", data = "230x80",
                              scale.data = "20x80", HVGs = "20"), class = "data.frame", row.names = "RNA"))
})

test_that("SeuratInfo works for version 5 object", {
  # checking that the last print() in SeuratInfo returns the correct value
  expect_equal(SeuratInfo(pbmc_small_v5),
               structure(list(default = "YES", counts = "230x80", data = "230x80",
                              scale.data = "20x80", HVGs = "20"), class = "data.frame", row.names = "RNA"))
})

