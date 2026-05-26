# sct2

Helper functions for single cell transcriptomics analysis in R.

## Installation

```r
devtools::install_github("pcantalupo/sct2")
```

## Usage

```r
library(sct2)

# Summarize a Seurat object
SeuratInfo(pbmc_small)

# Summarize Signac ChromatinAssays
SignacInfo(multiome_small)

# Find the metadata column matching the active ident
FindIdentLabel(pbmc_small)

# Save metadata to TSV
SaveMetadata(pbmc_small, file = "metadata.tsv")
```

## Command-line scripts

Scripts are in `inst/scripts/` and can be copied to `~/bin/` with `install_to_bin.sh`.

**seurat_info.R** — print a summary of a Seurat object (RDS or QS2)

```bash
seurat_info.R --seuratrds object.rds
seurat_info.R --seuratrds object.qs2
```

**seurat_save_metadata.R** — save Seurat metadata to a TSV file

```bash
seurat_save_metadata.R --seuratrds object.rds --outfile metadata.tsv
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
