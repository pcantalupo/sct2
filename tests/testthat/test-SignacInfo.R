test_that("SignacInfo works", {
  info = SignacInfo(multiome_small)

  expect_equal(info$ATAC$Cells, 120)
  expect_equal(info$ATAC$Features, 100)
  expect_equal(info$ATAC$Width_Range, "418-1423")
  expect_equal(info$ATAC$Width_Median, 912.5)
  expect_equal(info$ATAC$Fragments$n_fragment_objects, 6)
  expect_equal(info$ATAC$Fragments$fragment_num_cells, rep(20,6))  # 20 cells each sample
})

test_that("SignacInfo displays fragment paths for all 6 samples", {
  info = SignacInfo(multiome_small)

  paths = file.path("cellranger_count", c("2", "3", "4", "5a", "5b", "6"),
                    "outs", "atac_fragments.tsv.gz")

  expect_equal(info$ATAC$Fragments$fragment_paths, paths)
})
