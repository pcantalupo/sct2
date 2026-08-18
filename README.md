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

# Save metadata to TSV, or to RDS/QS2 to preserve factor levels and column classes
SaveMetadata(pbmc_small, file = "metadata.tsv")
SaveMetadata(pbmc_small, file = "metadata.rds")

# Read and write Seurat objects; format is inferred from the extension
WriteSeurat(pbmc_small, "pbmc_small.qs2")
seurat = ReadSeurat("pbmc_small.qs2")
```

## Command-line scripts

Scripts are in `inst/scripts/`. Symlink them from a directory on your `PATH` (e.g. `~/bin/`) to run them by name.

All scripts read `.rds`/`.RDS` or `.qs2`, inferred from the file extension.

**seurat_info.R** — print a summary of a Seurat object

```bash
seurat_info.R --seurat object.qs2
seurat_info.R --seurat object.qs2 --metadata   # also show metadata structure
```

**seurat_save_metadata.R** — save Seurat metadata to a TSV, RDS or QS2 file

```bash
seurat_save_metadata.R --seurat object.qs2 --outfile metadata.tsv
seurat_save_metadata.R --seurat object.qs2 --outfile metadata.rds   # preserves factor levels and column classes
```

**seurat_downsample.R** — randomly downsample cells and write a new object

```bash
seurat_downsample.R --seurat object.qs2 --downsample 0.05   # keep 5% of cells
seurat_downsample.R --seurat object.qs2 --downsample 50000  # keep 50k cells
```

Default output adds a `_ds<tag>` to the input basename and writes to the current directory (`--outfile` to override). `--downsample` is overloaded: a value of 1 or less is a fraction of cells to keep, a value above 1 is an absolute cell count. A fraction tags as a percentage (`0.3` → `_ds30`) and a count tags in thousands with a trailing `k` (`30` → `_ds0.03k`, `50000` → `_ds50k`), so the two never collide. A count larger than the object is an error rather than a no-op, since the `_ds<tag>` is resolved before the object is read and would otherwise name a downsample that never happened. Sampling uses a fixed `--seed` (default 1946). Use `--update` for objects written under an older SeuratObject, and `--force` to overwrite an existing output.

**seurat_update_object.R** — run `UpdateSeuratObject()` and save

```bash
seurat_update_object.R --seurat object.qs2 --outfile object_updated.qs2
```

Without `--outfile`, the input is overwritten in place. Because format is inferred from the extension, `--outfile` can also convert between `.rds` and `.qs2`. `UpdateSeuratObject()` validates the object it returns, so a failed migration errors before anything is written.

**seurat_strip_scaledata.R** — drop `scale.data` from every assay and save

```bash
seurat_strip_scaledata.R --seurat object.qs2                              # overwrite in place
seurat_strip_scaledata.R --seurat object.qs2 --outfile object_slim.qs2    # write elsewhere
```

Without `--outfile`, the input is overwritten in place; `--force` overwrites an existing `--outfile`. Because format is inferred from the extension, `--outfile` can also convert between `.rds` and `.qs2`. Removes all `scale.data*` layers from v5 assays and clears the `@scale.data` slot of v3 assays.

**seurat_dimplot_celltype-cluster.R** — UMAP DimPlot colored by celltype, labeled by a combined `<celltype>_<cluster>` column so every cluster of a celltype shares that celltype's color while staying individually labeled

```bash
seurat_dimplot_celltype-cluster.R --seurat object.qs2 --celltype singleR_cluster_labels --cluster RNA_snn_res.0.8
seurat_dimplot_celltype-cluster.R --seurat object.qs2 --downsample 0.1   # plot 10% of cells (or --downsample 2000 to plot 2000 cells)
```

Writes a PNG to `<outdir>/plots/` (`--outdir` defaults to `.`). Override the reduction with `--reduction`, the filename with `--outputfile`, the plot size with `--height`/`--width`, and pass `--repel` to repel the cluster labels. Use `--downsample` (fraction or absolute number of cells >= 2) to randomly downsample the Seurat object with a fixed `--seed`; cell count is noted in the plot subtitle.

**seurat_dimplot_splitby-colorby.R** — split UMAP DimPlot with one panel per `--splitby` value, colored by `--colorby` (mapped to `group.by`). The colorby column is coerced to a factor so every panel shares one color scale and a single unified legend.

```bash
seurat_dimplot_splitby-colorby.R --seurat object.qs2 --splitby RNA_snn_res.0.8 --colorby orig.ident
```

Writes a PNG to `plots/` by default. Override the reduction with `--reduction`, the output path with `--outputfile`, the plot size with `--height`/`--width`, and toggle cluster labels with `--label`/`--label_size`/`--repel`.

The panel grid is not configurable. `DimPlot_scCustom` assembles the per-level plots with `patchwork::wrap_plots()` and no `ncol`, so the layout falls through to `ggplot2:::wrap_dims()`, which calls `grDevices::n2mfrow(n)` on the number of `--splitby` levels. That function is landscape-biased, so it does not always produce a square-ish grid. Predict the layout without loading the object:

```bash
Rscript -e 'n <- 3; cat(rev(grDevices::n2mfrow(n)), "\n")'   # -> 1 3  (1 row, 3 columns)
```

| `--splitby` levels | rows × columns |
|---|---|
| 2 | 1 × 2 |
| 3 | 1 × 3 |
| 4 | 2 × 2 |
| 5 | 2 × 3 |
| 6 | 2 × 3 |
| 7–9 | 3 × 3 |

**seurat_dimplot_colorby.R** — one UMAP DimPlot per `--colorby` metadata column (mapped to `group.by`)

```bash
seurat_dimplot_colorby.R --seurat object.qs2 --colorby RNA_snn_res.0.8
seurat_dimplot_colorby.R --seurat object.qs2 --colorby orig.ident,RNA_snn_res.0.8,singleR_cluster_labels
```

`--colorby` takes a comma-separated list, so the object is loaded once and one PNG is written per column to `<outdir>/plots/` (`--outdir` defaults to `.`), named `<REDUCTION>_colored_by_<colorby>.png`. All columns are validated before the first plot is drawn. Non-factor colorby columns are coerced to a factor with naturally sorted levels. Override the reduction with `--reduction`, the plot size with `--height`/`--width`, and toggle labels with `--label`/`--label_size`/`--repel`.

**seurat_dotplot.R** — DotPlot of the top up-regulated genes per cluster from a markers table

```bash
seurat_dotplot.R --seurat object.qs2 --markers markers.rds --idents RNA_snn_res.0.8 --n_top_genes 5
```

`--markers` is a `FindAllMarkers()`-style RDS (`cluster`, `gene`, `avg_log2FC` columns). Set the idents with `--idents`, the number of genes per cluster with `--n_top_genes`, and the output PNG with `--output_path`. Adjust the plot with `--title`, `--labelsize`, `--rotatelabels`, `--width`, and `--height`.

## Functions

| Function | Description |
|---|---|
| `SeuratInfo()` | Summarize a Seurat object (idents, metadata, assays, reductions, graphs); prints the report and returns it invisibly as a `seurat_info` object |
| `print.seurat_info()` | Print method for the `seurat_info` object returned by `SeuratInfo()` |
| `SignacInfo()` | Summarize Signac ChromatinAssays within a Seurat object |
| `FindIdentLabel()` | Find the metadata column name that matches the active ident |
| `SaveMetadata()` | Save Seurat metadata to a TSV, RDS or QS2 file |
| `FixFragmentPaths()` | Fix paths to ATAC fragment files in a Seurat object |
| `FixClusterFactorLevels()` | Relevel cluster factors into numerical order |
| `DownsampleObject()` | Randomly subset cells from a Seurat or SingleCellExperiment object |
| `Self_scmapCluster()` | Map cell types within a single dataset using scmap |
| `TwoSample_scmapCluster()` | Map cell types between two datasets using scmap |
| `ReadSeurat()` | Read a Seurat object from an RDS or QS2 file |
| `WriteSeurat()` | Write a Seurat object to an RDS or QS2 file |
| `ValidateMetadataCols()` | Stop if any of the given metadata columns are missing from a Seurat object |
| `ValidateReduction()` | Stop, listing the available reductions, if a reduction is missing from a Seurat object |
