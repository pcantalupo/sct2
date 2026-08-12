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
annots_small = sample(annots, size = 100)
multiome_small@assays$ATAC@annotation = annots_small

# Keep 100 genes and 100 peaks to further reduce size of Seurat object
set.seed(42)
features_keep = c(sample(rownames(multiome_small), size = 100),
                  sample(rownames(multiome_small[['ATAC']]), size=100))
multiome_small = subset(multiome_small, features = features_keep)

# Replace the absolute fragment paths from the source object with neutral
# relative paths so no local directory structure ships with the package
frags = multiome_small@assays$ATAC@fragments
for (i in seq_along(frags)) {
  parts = strsplit(frags[[i]]@path, "[\\\\/]")[[1]]
  sample = parts[length(parts) - 2]
  slot(frags[[i]], "path") = file.path("cellranger_count", sample, "outs",
                                       "atac_fragments.tsv.gz")
}
multiome_small@assays$ATAC@fragments = frags

saveRDS(multiome_small, "inst/extdata/multiome_small.rds")
usethis::use_data(multiome_small, overwrite = TRUE)



