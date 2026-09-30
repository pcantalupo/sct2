# sct2

Helper functions for single cell transcriptomics analysis in R.

## Contents

- [Installation](#installation)
- [Usage](#usage)
- [Command-line scripts](#command-line-scripts)
  - [Object Info/Manipulation](#object-infomanipulation)
  - [Plotting](#plotting)
  - [Analysis](#analysis)
- [Functions](#functions)

## Installation

```r
devtools::install_github("pcantalupo/sct2")
```

## Usage

```r
library(sct2)

# Summarize a Seurat object
SeuratInfo(multiome_small)
```

```
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
```

```r
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

There are two general types of scripts: object information/manipulation and plotting. Get help for each script by passing `--help`

### Object Info/Manipulation

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

Default output adds a `_ds<tag>` to the input basename and writes to the current directory (`--outfile` to override). `--downsample` is overloaded: a value of 1 or less is a fraction of cells to keep, a value above 1 is an absolute cell count. A fraction tags as a percentage (`0.3` → `_ds30`) and a count tags in thousands with a trailing `k` (`30` → `_ds0.03k`, `50000` → `_ds50k`), so the two never collide. A count larger than the object is an error rather than a no-op, since the `_ds<tag>` is resolved before the object is read and would otherwise name a downsample that never happened.

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

### Plotting

All four plotting scripts share these options (pass the `--help` for more info):

| Option | Default | Notes |
|---|---|---|
| `--seurat` | required | Input object, `.rds` or `.qs2` |
| `--outdir` | `plots` | Written to directly; created after argument validation |
| `--outputfile` | derived | Filename only — a value containing `/` is rejected. Not available on `seurat_dimplot_colorby.R`, which writes one file per column |
| `--reduction` | `umap` | Not available on `seurat_dotplot.R` |
| `--width` | `7` | Inches |
| `--height` | `7` | Inches; `9` for `seurat_dotplot.R` |

Derived filenames start with `toupper(reduction)`, so `--reduction tsne` writes `TSNE_*`. Each script's derived name is listed with it below.

**seurat_dimplot_colorby.R** — one UMAP DimPlot per `--colorby` metadata column (mapped to `group.by`)

```bash
seurat_dimplot_colorby.R --seurat object.qs2 --colorby RNA_snn_res.0.8
seurat_dimplot_colorby.R --seurat object.qs2 --colorby orig.ident,RNA_snn_res.0.8,singleR_cluster_labels
```

`--colorby` takes a comma-separated list, so the object is loaded once and one PNG is written per column. Filename: `<REDUCTION>_colored_by_<colorby>.png`. All columns are validated before the first plot is drawn. Non-factor colorby columns are coerced to a factor with naturally sorted levels.

**seurat_dimplot_celltype-cluster.R** — UMAP DimPlot colored by celltype, labeled by a combined `<celltype>_<cluster>` column so every cluster of a celltype shares that celltype's color while staying individually labeled

```bash
seurat_dimplot_celltype-cluster.R --seurat object.qs2 --celltype singleR_cluster_labels --cluster RNA_snn_res.0.8
seurat_dimplot_celltype-cluster.R --seurat object.qs2 --downsample 0.1   # plot 10% of cells (or --downsample 2000 to plot 2000 cells)
```

Filename: `<REDUCTION>_colored_by_<celltype>_<cluster>.png`. Cell count is noted in the plot subtitle.

**seurat_dimplot_splitby-colorby.R** — split UMAP DimPlot with one panel per `--splitby` value, colored by `--colorby` (mapped to `group.by`). The colorby column is coerced to a factor so every panel shares one color scale and a single unified legend. Filename: `<REDUCTION>_split-<splitby>_color-<colorby>.png`.

```bash
seurat_dimplot_splitby-colorby.R --seurat object.qs2 --splitby RNA_snn_res.0.8 --colorby orig.ident
```

The panel grid is not configurable. `scCustomize::DimPlot_scCustom` assembles the per-level plots with `patchwork::wrap_plots()` and no `ncol`, so the layout falls through to `ggplot2:::wrap_dims()`. So it does not always produce a square-ish grid. The below table shows the grid layout for the number of `--splitby` levels

| `--splitby` levels | rows × columns |
|---|---|
| 2 | 1 × 2 |
| 3 | 1 × 3 |
| 4 | 2 × 2 |
| 5 | 2 × 3 |
| 6 | 2 × 3 |
| 7–9 | 3 × 3 |

**seurat_dotplot.R** — DotPlot of the top up-regulated genes per cluster from a markers table

```bash
seurat_dotplot.R --seurat object.qs2 --markers markers.rds --idents RNA_snn_res.0.8 --n_top_genes 5
```

`--markers` is a `FindAllMarkers()`-style RDS (`cluster`, `gene`, `avg_log2FC` columns). Set the idents with `--idents` and the number of genes per cluster with `--n_top_genes`. Filename: `dotplot_top<n_top_genes>_<idents>.png`.

### Analysis

**seurat_composition.R** — test how the composition of one metadata column (`--comp_var`, e.g. celltype) changes across the levels of another (`--group_var`, e.g. tissue_type), using samples as replicates

```bash
seurat_composition.R --input metadata.qs2 --group_var tissue_type --levels Tumor,Interface,Lung
seurat_composition.R --input object.qs2 --group_var tissue_type --sample_var patient --min_cells 50
```

`--input` is a Seurat object, a metadata data.frame from `seurat_save_metadata.R` (both `.rds` or `.qs2`), or a metadata `.tsv`/`.tsv.gz`. The metadata file loads much faster than the full object. Cells with NA in any of the three columns are dropped. Samples and comp levels keep their factor order; otherwise they are naturally sorted (`P2` before `P10`).

| Option | Default | Notes |
|---|---|---|
| `--group_var` | required | Column whose levels are compared |
| `--comp_var` | `celltype` | Column whose composition is measured |
| `--sample_var` | `patient` | Unit of replication |
| `--levels` | factor levels, else natural sort | Comma-separated order of `--group_var` levels; must list exactly the levels in the data |
| `--min_cells` | `0` | Drop sample x group combinations with fewer cells (logged). Dropped cells are left out of every output |
| `--transform` | `logit` | `logit` or `asin` (arcsine square root) |
| `--outdir` | `plots` | PNG directory |
| `--resultsdir` | `results` | TSV and XLSX directory |

Statistics follow the propeller method from the Bioconductor package `speckle`: each sample's proportions are transformed, then limma tests every pair of `--group_var` levels, earlier level in `--levels` as the numerator. When a sample appears in more than one group, the fit blocks on sample with `limma::duplicateCorrelation()`, as the speckle vignette recommends for repeated measures. The log says which fit was used, and warns when a group has only one sample. Pooled cell counts are never tested. The values limma fits are saved in the `transformed_prop` column of the `_per_sample` table.

Output files are named `composition_<group_var>_<comp_var>_*`; characters other than `A-Za-z0-9._-` become `_`.

| Output | Contents |
|---|---|
| `_pooled.tsv/.xlsx` | Cells pooled across samples. One row per comp level plus a `Total` row. For each group: `<group>_n`, `<group>_pct_of_group` (share of the group's cells) and `<group>_pct_of_comp` (share of the comp level's cells). `Total_n` is the comp level's cell count |
| `_per_sample.tsv/.xlsx` | One row per sample x group x comp level: `n`, `total`, `prop`, and `transformed_prop`, the value limma fits. Zero counts are included. Logit adds 0.5 to each count first, so `transformed_prop` is not exactly `logit(prop)` |
| `_enrichment.tsv/.xlsx` | Pooled `observed`, `expected` (row total x column total / grand total) and `log2_enrichment`, which is log2 of Ro/e. Descriptive only |
| `_stats.tsv/.xlsx` | One row per contrast x comp level: mean per-sample proportion in each group, `prop_ratio` (their ratio; `Inf` or `0` when one mean is 0), `estimate` (difference on the transformed scale; natural-log odds ratio for logit), `t`, `p_value`, `fdr` (BH within each contrast), `model`, `transform` |
| `_pooled_stackedbar.png` | Pooled composition per group |
| `_per_sample_stackedbar.png` | Composition per sample, faceted by group |
| `_per_sample_stackedbar_by_sample.png` | Same bars faceted by sample, so a sample's groups sit side by side |
| `_per_sample_boxplot.png` | Per-sample proportions by group, one panel per comp level; lines join a sample's points |
| `_enrichment_heatmap.png` | `log2_enrichment`; grey where a comp level has no cells in a group |

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
