library(Seurat)

dir.create("inst/extdata", showWarnings = FALSE)
saveRDS(pbmc_small, "inst/extdata/pbmc_small.rds")
usethis::use_data(pbmc_small, overwrite = TRUE)

