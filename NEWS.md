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
