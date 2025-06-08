test_that("SignacInfo works", {
  info = capture.output(SignacInfo(multiome_small))
  expect_equal(info[4],
               "ATAC: 53984 peaks, 120 cells, width range [157-3228], median 929, fragments: 6")
})

test_that("SignacInfo displays fragment paths for all 6 samples", {
  info = capture.output(SignacInfo(multiome_small))

  paths = c("  C:\\Users\\pgc92\\OneDrive - University of Pittsburgh\\DBMI_Genomics_Analysis_Chandran\\Paul\\dasilva_multiome\\cellranger_count\\2\\outs\\atac_fragments.tsv.gz (20 cells)",
            "  C:\\Users\\pgc92\\OneDrive - University of Pittsburgh\\DBMI_Genomics_Analysis_Chandran\\Paul\\dasilva_multiome\\cellranger_count\\3\\outs\\atac_fragments.tsv.gz (20 cells)",
            "  C:\\Users\\pgc92\\OneDrive - University of Pittsburgh\\DBMI_Genomics_Analysis_Chandran\\Paul\\dasilva_multiome\\cellranger_count\\4\\outs\\atac_fragments.tsv.gz (20 cells)",
            "  C:\\Users\\pgc92\\OneDrive - University of Pittsburgh\\DBMI_Genomics_Analysis_Chandran\\Paul\\dasilva_multiome\\cellranger_count\\5a\\outs\\atac_fragments.tsv.gz (20 cells)",
            "  C:\\Users\\pgc92\\OneDrive - University of Pittsburgh\\DBMI_Genomics_Analysis_Chandran\\Paul\\dasilva_multiome\\cellranger_count\\5b\\outs\\atac_fragments.tsv.gz (20 cells)",
            "  C:\\Users\\pgc92\\OneDrive - University of Pittsburgh\\DBMI_Genomics_Analysis_Chandran\\Paul\\dasilva_multiome\\cellranger_count\\6\\outs\\atac_fragments.tsv.gz (20 cells)")

  expect_equal(info[5:10], paths)
})
