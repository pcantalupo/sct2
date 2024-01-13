test_that("TwoSample_scmapCluster works", {
  twosample = TwoSample_scmapCluster(pbmc_small_sce, pbmc_small_sce)

  versusnames = c("1vs2", "2vs1")
  expect_equal(names(twosample), versusnames)

  clustercol = "RNA_snn_res.0.8"
  expect_equal(names(twosample[[versusnames[1]]]), clustercol)

  listelements = c("scmap_cluster_labs", "scmap_cluster_siml", "combined_labs")
  expect_equal(names(twosample[[versusnames[1]]][[clustercol]]),
               listelements)

  expect_equal(names(twosample[[versusnames[2]]]), clustercol)
})
