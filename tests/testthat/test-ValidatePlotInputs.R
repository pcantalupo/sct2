test_that("ValidateMetadataCols passes silently when all columns exist", {
  expect_null(ValidateMetadataCols(pbmc_small, c("orig.ident", "groups")))
})

test_that("ValidateMetadataCols stops and lists every missing column", {
  expect_error(ValidateMetadataCols(pbmc_small, c("orig.ident", "nope1", "nope2")),
               "nope1, nope2")
})

test_that("ValidateReduction passes silently when the reduction exists", {
  expect_null(ValidateReduction(pbmc_small, "pca"))
})

test_that("ValidateReduction stops and lists the available reductions", {
  # pbmc_small_v5 ships only the tsne reduction, no umap
  expect_error(ValidateReduction(pbmc_small_v5, "umap"), "tsne")
})
