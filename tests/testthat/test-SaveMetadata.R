test_that("SaveMetadata writes a metadata file", {
  file = tempfile(fileext = ".tsv")
  on.exit(unlink(file))
  SaveMetadata(seurat = pbmc_small, filename = file)
  #Test if the file was created
  expect_true(file.exists(file))
})

test_that("SaveMetadata uses the colname_for_rows parameter appropriately", {
  file = tempfile(fileext = ".tsv")
  on.exit(unlink(file))
  colname = "foo"
  SaveMetadata(seurat = pbmc_small, filename = file, colname_for_rows = colname)
  metadata = read.delim(file, sep="\t")
  expect_equal(colnames(metadata)[1], colname)
})

test_that("SaveMetadata writes an RDS file that preserves factors and rownames", {
  file = tempfile(fileext = ".rds")
  on.exit(unlink(file))
  seurat = pbmc_small
  seurat$mycluster = factor(seurat$RNA_snn_res.1, levels = c("2", "0", "1"))
  SaveMetadata(seurat = seurat, filename = file)

  metadata = readRDS(file)
  expect_s3_class(metadata$mycluster, "factor")
  expect_equal(levels(metadata$mycluster), c("2", "0", "1"))
  expect_equal(rownames(metadata), colnames(seurat))
  expect_equal(metadata, seurat[[]])
})

test_that("SaveMetadata writes a QS2 file that preserves factors and rownames", {
  skip_if_not_installed("qs2")
  file = tempfile(fileext = ".qs2")
  on.exit(unlink(file))
  seurat = pbmc_small
  seurat$mycluster = factor(seurat$RNA_snn_res.1, levels = c("2", "0", "1"))
  SaveMetadata(seurat = seurat, filename = file)

  metadata = qs2::qs_read(file)
  expect_s3_class(metadata$mycluster, "factor")
  expect_equal(levels(metadata$mycluster), c("2", "0", "1"))
  expect_equal(rownames(metadata), colnames(seurat))
  expect_equal(metadata, seurat[[]])
})

test_that("SaveMetadata expects seurat and filename params", {
  expect_error(SaveMetadata(), "seurat and filename parameters are required")
  expect_error(SaveMetadata(seurat = pbmc_small), "seurat and filename parameters are required")
  expect_error(SaveMetadata(filename = "pbmc_small_metadata.tsv"), "seurat and filename parameters are required")
})

