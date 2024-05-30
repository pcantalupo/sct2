library(Seurat)

dir.create("inst/extdata", showWarnings = FALSE)

# Seurat Version 4
#
saveRDS(pbmc_small, "inst/extdata/pbmc_small.rds")
usethis::use_data(pbmc_small, overwrite = TRUE)   # saves the .rda file by same name as the object

pbmc_small_sce = as.SingleCellExperiment(pbmc_small)
saveRDS(pbmc_small_sce, "inst/extdata/pbmc_small_sce.rds")
usethis::use_data(pbmc_small_sce, overwrite = TRUE)


# Seurat Version 5
#
counts = GetAssayData(pbmc_small, layer = "counts")
pbmc_small_v5 = CreateSeuratObject(counts, meta.data = pbmc_small[[]])
LayerData(pbmc_small_v5, layer = "data") = GetAssayData(pbmc_small, layer = "data")
LayerData(pbmc_small_v5, layer = "scale.data") = GetAssayData(pbmc_small, layer = "scale.data")
VariableFeatures(pbmc_small_v5) = VariableFeatures(pbmc_small)
pbmc_small_v5[['tsne']] = pbmc_small[['tsne']]
saveRDS(pbmc_small_v5, "inst/extdata/pbmc_small_v5.rds")
usethis::use_data(pbmc_small_v5, overwrite = TRUE)

