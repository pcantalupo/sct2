library(Seurat)
library(Signac)

dir.create("inst/extdata", showWarnings = FALSE)

multiome = readRDS("~/Desktop/seurat_harmony.rds")
Idents(multiome) = "orig.ident"

multiome_small = subset(multiome, downsample = 20)
saveRDS(multiome_small, "inst/extdata/multiome_small.rds")
usethis::use_data(multiome_small, overwrite = TRUE)



