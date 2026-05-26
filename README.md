# sct2

Helper functions for single cell transcriptomics analysis in R.

## Installation

```r
devtools::install_github("pcantalupo/sct2")
```

## Functions

| Function | Description |
|---|---|
| `SeuratInfo()` | Summarize a Seurat object (idents, metadata, assays, reductions, graphs) |
| `SignacInfo()` | Summarize Signac ChromatinAssays within a Seurat object |
| `FindIdentLabel()` | Find the metadata column name that matches the active ident |
| `SaveMetadata()` | Save Seurat metadata to a TSV file |
| `FixFragmentPaths()` | Fix paths to ATAC fragment files in a Seurat object |
| `FixClusterFactorLevels()` | Relevel cluster factors into numerical order |
| `SelfScmapCluster()` | Map cell types within a single dataset using scmap |
| `TwoSampleScmapCluster()` | Map cell types between two datasets using scmap |
