
seurat = readRDS("atomx_small.rds")
# RNA assay has data slot BUT it is NOT normalized
#     LayerData(seurat, layer = "data", assay = "RNA")[1:10,1:5]
# MettoLungBoneRNA assay contains the normalized data
#     LayerData(seurat, layer = "data", assay = "Met.to.Lung.Bone_Normalization.RNA.1_1")[1:10,1:5]
Seurat::DefaultAssay(seurat) = "Met.to.Lung.Bone_Normalization.RNA.1_1"

niches = "Met.to.Lung.Bone_Neighborhood.Analysis.1_1_assignments"
celltypes = "RNA_Met.to.Lung.Bone_Cell.Typing.InSituType.1_1_clusters"
clusters = "Met.to.Lung.Bone_Neighbor.network.expression.space.1_1_cluster_Met.to.Lung.Bone_Leiden.Clustering.1_1"
Seurat::Idents(seurat) = celltypes

start <- Sys.time()
summarize_seurat(seurat)
end <- Sys.time()
end - start

start <- Sys.time()
SeuratInfo(seurat)
end <- Sys.time()
end - start



seurat

