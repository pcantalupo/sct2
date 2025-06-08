library(Seurat)
library(Signac)

dir.create("inst/extdata", showWarnings = FALSE)

# Load Multiome Seurat object
multiome = readRDS("~/Desktop/seurat_harmony.rds")
Idents(multiome) = "orig.ident"

# Downsample to 20 cells per sample (120 cells total)
multiome_small = subset(multiome, downsample = 20)

# Reduce size of the Annotation from 118 MB to 42 KB by keeping just 100 annotation rows
annots = Annotation(multiome_small[['ATAC']])
annots_small = sample(granges(annots), size = 100)
multiome_small@assays$ATAC@annotation = annots_small

# Keep 100 genes and 100 peaks to further reduce size of Seurat object
set.seed(42)
features_keep = c(sample(rownames(multiome_small), size = 100),
                  sample(rownames(multiome_small[['ATAC']]), size=100))
multiome_small = subset(multiome_small, features = features_keep)

saveRDS(multiome_small, "inst/extdata/multiome_small.rds")
usethis::use_data(multiome_small, overwrite = TRUE)



