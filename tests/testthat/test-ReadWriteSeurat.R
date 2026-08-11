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
  outdir <- tempfile()
  dir.create(outdir)
  on.exit(unlink(outdir, recursive = TRUE))
  path <- file.path(outdir, "seurat.rds")
  WriteSeurat(pbmc_small, path)
  expect_equal(list.files(outdir), basename(path))
})

# The atomic write only shows up when serialization fails partway through: on
# success the tempfile has already been renamed away, so a direct saveRDS() to
# path is indistinguishable. These mock the writer to reproduce that failure --
# it must clobber whatever file it was handed *before* throwing, the way a real
# truncating write does, or the mock passes against a non-atomic implementation
# that never touched the destination.

fail_after_clobbering = function(object, file, ...) {
  writeBin(raw(1), file)
  stop("boom")
}

test_that("WriteSeurat leaves the existing rds intact when serialization fails", {
  outdir = tempfile()
  dir.create(outdir)
  on.exit(unlink(outdir, recursive = TRUE))
  path = file.path(outdir, "seurat.rds")
  WriteSeurat(pbmc_small, path)
  before = unname(tools::md5sum(path))

  local_mocked_bindings(saveRDS = fail_after_clobbering, .package = "base")
  expect_error(WriteSeurat(pbmc_small, path), "boom")
  expect_equal(unname(tools::md5sum(path)), before)
  expect_equal(list.files(outdir), basename(path))
})

test_that("WriteSeurat leaves the existing qs2 intact when serialization fails", {
  outdir = tempfile()
  dir.create(outdir)
  on.exit(unlink(outdir, recursive = TRUE))
  path = file.path(outdir, "seurat.qs2")
  WriteSeurat(pbmc_small, path)
  before = unname(tools::md5sum(path))

  local_mocked_bindings(qs_save = fail_after_clobbering, .package = "qs2")
  expect_error(WriteSeurat(pbmc_small, path), "boom")
  expect_equal(unname(tools::md5sum(path)), before)
  expect_equal(list.files(outdir), basename(path))
})
