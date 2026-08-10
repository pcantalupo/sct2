# sct2 0.3.7

- `seurat_downsample.R` writes its default `_ds<tag>` output to the current directory instead of next to the input object
- Document the panel grid layout of `seurat_dimplot_splitby-colorby.R`
- Add Claude Code GitHub Workflow

# sct2 0.3.6

- Add `seurat_dimplot_colorby.R` command-line script: one DimPlot per `--colorby` metadata column, from a single object load
- Add `CLAUDE.md` documenting the bundled test objects and the command-line script conventions

# sct2 0.3.5

- `SaveMetadata()` picks the output format from the filename extension: `.rds` and `.qs2` save the metadata data.frame directly, preserving factor levels, column classes, and cell barcodes as rownames; any other extension writes a TSV as before. `colname_for_rows` applies to TSV output only
- Add `qs2` and `tools` to DESCRIPTION Imports

# sct2 0.3.4

- Add `seurat_strip_scaledata.R` command-line script
- Add `--downsample`, `--ncells`, and `--seed` flags to `seurat_dimplot_celltype-cluster.R` to plot a random subset of cells without subsetting the object
- Add a test covering `FixClusterFactorLevels` releveling of out-of-order cluster levels

# sct2 0.3.3

- Add `seurat_dimplot_splitby-colorby.R` command-line script: split UMAP DimPlot (one panel per `--splitby` value) colored by a `--colorby` metadata column
- Add progress messages to `seurat_dotplot.R`
- Remove `install_to_bin.sh`; scripts are now symlinked from `bin/` directly

# sct2 0.3.2

- Add `seurat_dotplot.R` command-line script: DotPlot of the top up-regulated genes per cluster from a markers table

# sct2 0.3.1

- Add `seurat_dimplot_celltype-cluster.R` command-line script: UMAP DimPlot colored by celltype, labeled by a combined `<celltype>_<cluster>` column

# sct2 0.3.0

- Add `seurat_downsample.R` and `seurat_update_object.R` command-line scripts
- Add `ReadSeurat()` and `WriteSeurat()` functions and refactor scripts to use them
- `ReadSeurat()` now errors if the loaded object is not a Seurat object
- Standardize command-line script flags to `--seurat` (input) and `--outfile` (output) across all four scripts, replacing `--seuratrds`, `--input`, and `--output`

# sct2 0.2.3

- Fix `FixClusterFactorLevels` to skip non-factor `snn_res` columns (e.g. score columns)

# sct2 0.2.2

- Add `--metadata` flag to `seurat_info.R` to display metadata via `SeuratInfo()`
- Add qs2 support to `seurat_save_metadata.R`
- Add `.Rbuildignore` and `NEWS.md`
- Update README with usage examples and command-line scripts

# sct2 0.2.1

- Fix `FindIdentLabel` error when `seurat@meta.data` contains POSIXct columns (#1)

# sct2 0.2.0

- Add qs2 support
- Implement faster `SeuratInfo` algorithm
- Complete rewrite of `SignacInfo`

# sct2 0.1.0

- Add `FixFragmentPaths` function
- Add `FixClusterFactorLevels` function
- Add `SaveMetadata` function
- Add `SelfScmapCluster` and `TwoSampleScmapCluster` functions
- Add `FindIdentLabel` function
