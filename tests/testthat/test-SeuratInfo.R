test_that("SeuratInfo reports assays, reductions and idents for a v3 object", {
  expect_snapshot(SeuratInfo(pbmc_small))
})

test_that("SeuratInfo reports the same summary for a v5 object", {
  expect_snapshot(SeuratInfo(pbmc_small_v5))
})

test_that("SeuratInfo shows the metadata when asked", {
  expect_snapshot(SeuratInfo(pbmc_small, metadata = TRUE))
})

test_that("SeuratInfo reports both assays of a multiome object", {
  expect_snapshot(SeuratInfo(multiome_small))
})

test_that("SeuratInfo returns a seurat_info object invisibly", {
  expect_invisible(SeuratInfo(pbmc_small))

  info = SeuratInfo(pbmc_small)
  expect_s3_class(info, "seurat_info")
  expect_equal(info$version, "4.0.0")
  expect_null(info$metadata)
  expect_equal(info$graphs, "RNA_snn")
  expect_equal(info$reductions, c(pca = "RNA", tsne = "RNA"))
  expect_equal(info$ident_label, "RNA_snn_res.1")
  expect_equal(as.integer(info$idents), c(36L, 25L, 19L))
})

test_that("SeuratInfo preserves the 0.4.6 assay data.frame as the assays element", {
  info = SeuratInfo(pbmc_small)
  expect_identical(info$assays,
                   structure(list(default = "YES",
                                  counts = "230x80",
                                  data = "230x80",
                                  scale.data = "20x80",
                                  HVGs = "20"),
                             class = "data.frame",
                             row.names = "RNA"))
})

test_that("SeuratInfo stores the metadata when metadata = TRUE", {
  info = SeuratInfo(pbmc_small, metadata = TRUE)
  expect_identical(info$metadata, pbmc_small[[]])
})

test_that("SeuratInfo reports the reductions and assays of a multiome object", {
  info = SeuratInfo(multiome_small)
  expect_equal(info$reductions, c(SCT_pca_umap = "SCT", umap = "SCT"))
  expect_equal(rownames(info$assays), c("RNA", "ATAC"))
  expect_equal(info$assays$default, c("YES", ""))
})
