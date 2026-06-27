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

All scripts read `.rds`/`.RDS` or `.qs2`, inferred from the file extension.

**seurat_info.R** — print a summary of a Seurat object

```bash
seurat_info.R --seurat object.qs2
seurat_info.R --seurat object.qs2 --metadata   # also show metadata structure
```

**seurat_save_metadata.R** — save Seurat metadata to a TSV file

```bash
seurat_save_metadata.R --seurat object.qs2 --outfile metadata.tsv
```

**seurat_downsample.R** — randomly downsample cells and write a new object

```bash
seurat_downsample.R --seurat object.qs2 --downsample 0.05   # keep 5% of cells
seurat_downsample.R --seurat object.qs2 --ncells 50000      # keep 50k cells
```

Default output adds a `_ds<tag>` next to the input (`--outfile` to override). Use `--update` for objects written under an older SeuratObject, and `--force` to overwrite an existing output.

**seurat_update_object.R** — run `UpdateSeuratObject()` and save

```bash
seurat_update_object.R --seurat object.qs2 --outfile object_updated.qs2
```

Without `--outfile`, the input is overwritten in place. Because format is inferred from the extension, `--outfile` can also convert between `.rds` and `.qs2`.

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
