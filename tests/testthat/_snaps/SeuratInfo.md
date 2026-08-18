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

