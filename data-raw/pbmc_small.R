library(Seurat)

dir.create("inst/extdata", showWarnings = FALSE)
saveRDS(pbmc_small, "inst/extdata/pbmc_small.rds")
usethis::use_data(pbmc_small, overwrite = TRUE)

pbmc_small_sce = as.SingleCellExperiment(pbmc_small)
saveRDS(pbmc_small_sce, "inst/extdata/pbmc_small_sce.rds")
usethis::use_data(pbmc_small_sce, overwrite = TRUE)

