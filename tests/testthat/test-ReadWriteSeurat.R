test_that("WriteSeurat and ReadSeurat round-trip via rds", {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  WriteSeurat(pbmc_small, path)
  result <- ReadSeurat(path)
  expect_s4_class(result, "Seurat")
  expect_equal(ncol(result), ncol(pbmc_small))
})

test_that("WriteSeurat and ReadSeurat round-trip via qs2", {
  path <- tempfile(fileext = ".qs2")
  on.exit(unlink(path))
  WriteSeurat(pbmc_small, path)
  result <- ReadSeurat(path)
  expect_s4_class(result, "Seurat")
  expect_equal(ncol(result), ncol(pbmc_small))
})

test_that("ReadSeurat errors on unsupported extension", {
  expect_error(ReadSeurat("object.txt"), "Unsupported")
})

test_that("ReadSeurat errors when the file is not a Seurat object", {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(data.frame(x = 1:3), path)
  expect_error(ReadSeurat(path), "not a Seurat object")
})

test_that("WriteSeurat errors on unsupported extension", {
  expect_error(WriteSeurat(pbmc_small, "object.txt"), "Unsupported")
})

test_that("WriteSeurat can overwrite its own input path in place", {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  WriteSeurat(pbmc_small, path)

  seurat <- ReadSeurat(path)
  WriteSeurat(seurat, path)

  result <- ReadSeurat(path)
  expect_s4_class(result, "Seurat")
  expect_equal(ncol(result), ncol(pbmc_small))
})

test_that("WriteSeurat leaves no stray tempfile in the target directory", {
  dir <- tempfile()
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  path <- file.path(dir, "seurat.rds")
  WriteSeurat(pbmc_small, path)
  expect_equal(list.files(dir), basename(path))
})
