test_that("SaveMetadata writes a metadata file", {
  file = "pbmc_small_metadata.tsv"
  SaveMetadata(seurat = pbmc_small, filename = file)
  #Test if the file was created
  expect_true(file.exists(file))
  unlink(file)
})

test_that("SaveMetadata uses the colname_for_rows parameter appropriately", {
  file = "pbmc_small_metadata.tsv"
  colname = "foo"
  SaveMetadata(seurat = pbmc_small, filename = file, colname_for_rows = colname)
  metadata = read.delim(file, sep="\t")
  expect_equal(colnames(metadata)[1], colname)
  unlink(file)
})

test_that("SaveMetadata expects seurat and filename params", {
  expect_error(SaveMetadata(), "seurat and filename parameters are required")
  expect_error(SaveMetadata(seurat = pbmc_small), "seurat and filename parameters are required")
  expect_error(SaveMetadata(filename = "pbmc_small_metadata.tsv"), "seurat and filename parameters are required")
})

