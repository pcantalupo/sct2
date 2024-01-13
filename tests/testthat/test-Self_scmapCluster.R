test_that("Self_scmapCluster works", {
  selfscmap = Self_scmapCluster(pbmc_small_sce)

  clustercol = "RNA_snn_res.0.8"
  expect_equal(names(selfscmap), clustercol)

  listelements = c("scmap_cluster_labs", "scmap_cluster_siml", "combined_labs")
  expect_equal(names(selfscmap[[clustercol]]), listelements)
})
