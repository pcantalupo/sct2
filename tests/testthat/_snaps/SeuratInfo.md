# SeuratInfo reports assays, reductions and idents for a v3 object

    Code
      SeuratInfo(pbmc_small)
    Output
      Seurat version:  4.0.0 
      
      Graphs: RNA_snn
      Reductions: pca (RNA), tsne (RNA)
      Ident label: RNA_snn_res.1
      Idents():
             0  1  2
      Count 36 25 19
      
      Assays:
          default counts   data scale.data HVGs
      RNA     YES 230x80 230x80      20x80   20

# SeuratInfo reports the same summary for a v5 object

    Code
      SeuratInfo(pbmc_small_v5)
    Output
      Seurat version:  5.0.2 
      
      Graphs: 
      Reductions: tsne (RNA)
      Ident label: orig.ident
      Idents():
            SeuratProject
      Count            80
      
      Assays:
          default counts   data scale.data HVGs
      RNA     YES 230x80 230x80      20x80   20

# SeuratInfo shows the metadata when asked

    Code
      SeuratInfo(pbmc_small, metadata = TRUE)
    Output
      Seurat version:  4.0.0 
      
      Metadata: 'data.frame':	80 obs. of  7 variables:
       $ orig.ident     : Factor w/ 1 level "SeuratProject": 1 1 1 1 1 1 1 1 1 1 ...
       $ nCount_RNA     : num  70 85 87 127 173 70 64 72 52 100 ...
       $ nFeature_RNA   : int  47 52 50 56 53 48 36 45 36 41 ...
       $ RNA_snn_res.0.8: Factor w/ 2 levels "0","1": 1 1 2 1 1 1 1 1 1 1 ...
       $ letter.idents  : Factor w/ 2 levels "A","B": 1 1 2 1 1 1 1 1 1 1 ...
       $ groups         : chr  "g2" "g1" "g2" "g2" ...
       $ RNA_snn_res.1  : Factor w/ 3 levels "0","1","2": 1 1 1 1 1 1 1 1 1 1 ...
      
      Graphs: RNA_snn
      Reductions: pca (RNA), tsne (RNA)
      Ident label: RNA_snn_res.1
      Idents():
             0  1  2
      Count 36 25 19
      
      Assays:
          default counts   data scale.data HVGs
      RNA     YES 230x80 230x80      20x80   20

# SeuratInfo reports both assays of a multiome object

    Code
      SeuratInfo(multiome_small)
    Output
      Seurat version:  5.1.0 
      
      Graphs: 
      Reductions: SCT_pca_umap (SCT), umap (SCT)
      Ident label: orig.ident
      Idents():
             2  3  4 5a 5b  6
      Count 20 20 20 20 20 20
      
      Assays:
           default  counts    data scale.data HVGs
      RNA      YES 100x120 100x120      7x120    7
      ATAC         100x120 100x120        0x0    0

