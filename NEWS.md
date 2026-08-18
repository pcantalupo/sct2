# sct2 0.4.7

- Split `SeuratInfo()` into data gathering and formatting. It now builds a `seurat_info` classed list holding the version, metadata, graphs, reductions, images, ident label, idents table and assay table, prints it, and returns it invisibly. All formatting lives in the new exported `print.seurat_info()`. The printed report is unchanged, so every call site behaves as before
- **Breaking**: `SeuratInfo()` returns a `seurat_info` list instead of the assay `data.frame`. The former return value is preserved verbatim as the `assays` element, so `SeuratInfo(x)` becomes `SeuratInfo(x)$assays`
- Add structure assertions to `test-SeuratInfo.R`, plus snapshots for the `metadata = TRUE` branch and for `multiome_small`, which covers the multi-assay table and reductions reporting an assay absent from the object

# sct2 0.4.6

- `SeuratInfo()` ends with `invisible(df)` and documents a `@return` tag, so the assay table is a declared return value rather than a side effect of the final `print.data.frame` call (#7)
- Convert `test-SeuratInfo.R` to `expect_snapshot()`. The two `expect_equal()` assertions covered only the assay table; the snapshot captures the whole printed report, including the version line, graphs, reductions, ident label and `Idents()` counts, none of which were tested before

# sct2 0.4.5

- Add `ValidateMetadataCols()` and `ValidateReduction()`, and call them from all three `seurat_dimplot_*.R` scripts in place of each script's own copy of the check (#12). `seurat_dimplot_celltype-cluster.R` previously hand-rolled a per-column `%in%` test with its own message format and reported only the first missing column; it now reports them all at once in the shared format
- `seurat_dimplot_celltype-cluster.R` and `seurat_dimplot_splitby-colorby.R` validate `--reduction` against the object's reductions. Neither checked it before, so a bad value failed later inside `DimPlot()`

# sct2 0.4.4

- Add `LexSCT_snn_res.0.8` to `multiome_small`: the same cluster assignments as `SCT_snn_res.0.8`, but with its 18 factor levels in lexicographic order. Every other cluster column in the object already has numerically ordered levels, so there was no fixture for the case `FixClusterFactorLevels()` exists to fix
- Rewrite the first test in `test-FixClusterFactorLevels.R` to use `multiome_small` and assert the expected level vectors as literals (#8). It previously used `pbmc_small`, whose three clusters cannot distinguish lexicographic from numeric order, and compared the returned levels to their own sorted version, so it passed against a no-op implementation
- Document `multiome_small`'s assays, reductions, metadata columns and cluster level ordering

# sct2 0.4.3

- `seurat_dimplot_celltype-cluster.R` uses `DownsampleObject()` instead of its own copy of the sampling logic (#9). Both scripts that downsample now share the exported function
- Remove the `--ncells` flag from `seurat_dimplot_celltype-cluster.R`; pass an absolute count to `--downsample` instead. The default `--downsample` is 1 (no downsampling) and the default `--seed` is 1946
- `seurat_dimplot_celltype-cluster.R` subsets the object rather than passing `cells =` to `DimPlot()`. The celltype color mapping is built from the full celltype set before the downsample, so a celltype that loses all its cells no longer shifts every subsequent celltype's position in `colors_polychrome`

# sct2 0.4.2

- `seurat_downsample.R` uses `DownsampleObject()` instead of its own copy of the sampling logic (#9). `seurat_dimplot_celltype-cluster.R` still carries its own copy
- `seurat_downsample.R` forwards `--seed` to the sampler. It was parsed and printed but never passed, so every run sampled with the default seed regardless of the flag
- `seurat_downsample.R` stops when an absolute `--downsample` count exceeds the object size, rather than relying on the clamp in `DownsampleObject()`. The `_ds<tag>` output name is resolved before the object is read, so a silent clamp would write a full-size object under a name claiming a downsample that never happened
- `seurat_downsample.R` rejects a `--downsample` that is not positive, and `DownsampleObject()` errors when a fraction floors to zero cells instead of failing later inside `subset.Seurat()` with "No cells found"
- `seurat_downsample.R` output tags distinguish a fraction from a count: a fraction tags as a percentage and a count tags in thousands with a trailing `k`, so `--downsample 0.3` gives `_ds30` and `--downsample 30` gives `_ds0.03k`. Both produced `_ds30` before. `--downsample 1` tags `_ds100` rather than `_ds1`
- `seurat_downsample.R` logs the before and after cell counts with the seed used
- Remove the `--ncells` flag from `seurat_downsample.R`; pass an absolute count to `--downsample` instead. The default `--seed` is 1946
- Add `test-seurat_downsample.R` covering the output tag scheme and the argument guards

# sct2 0.4.1

- Add `DownsampleObject()`: draw a uniform random subset of cells from a Seurat or SingleCellExperiment object. `downsample` is a fraction when `<= 1` and an absolute cell count when `> 1`, clamped to the object size; the object is returned unchanged and the RNG untouched when the resolved count is the whole object (#9)
- Scope the `SaveMetadata` test file cleanup with `on.exit()`
- Replace personal fragment paths in `multiome_small` with relative paths
- Document that the scmap tests are out of scope for testing work

# sct2 0.4.0

- Remove the `--check` flag from `seurat_update_object.R`. `UpdateSeuratObject()` runs `validObject()` on the object and on every image before returning, so a failed migration already errors before anything is written (#6). Invocations passing `--check` will now fail with an invalid-flag error
- `seurat_update_object.R` validates `--outfile` in the options block, before the object load: the extension must be `.rds` or `.qs2`, and the parent directory must exist and be writable (#5). Previously a bad output path was not discovered until after the load and update
- `WriteSeurat()` writes atomically: the object is serialized to a tempfile in the target directory and renamed into place, so an interrupt or a failed write can no longer corrupt an existing file at that path (#4). This matters for the scripts that overwrite their input in place. Note that an in-place rewrite now needs free space for a second copy, and write permission on the directory rather than on the file
- Add tests covering the `--outfile` guards and the failed-write path, and require `testthat (>= 3.2.0)`

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
