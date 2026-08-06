# sct2

Helper functions and CLI scripts for single cell transcriptomics analysis. Exported functions live in `R/`, CLI scripts in `inst/scripts/`.

## Test data

The package ships Seurat objects for testing — use these instead of asking for a file:

| Object | Notes |
|---|---|
| `pbmc_small` | v3 Seurat, 230 genes x 80 cells. Reductions: `pca`, `tsne` |
| `pbmc_small_v5` | same object as a v5 Seurat. Reduction: `tsne` only — **no `umap`** |
| `pbmc_small_sce` | `pbmc_small` as a SingleCellExperiment |
| `multiome_small` | 6 samples x 20 cells, RNA + ATAC (Signac `ChromatinAssay`) |

Metadata columns on the pbmc objects: `orig.ident`, `nCount_RNA`, `nFeature_RNA`, `RNA_snn_res.0.8`, `RNA_snn_res.1`, `letter.idents`, `groups`.

## Testing a CLI script

Scripts call `pacman::p_load(sct2)`, which loads the *installed* package — `devtools::load_all()` does not reach them. Write the test object to a temp `.qs2` with `devtools::load_all()` in a separate `Rscript -e`, then invoke the script against that file:

```bash
Rscript -e "suppressMessages(devtools::load_all('.')); WriteSeurat(pbmc_small_v5, '/tmp/pbmc.qs2')"
inst/scripts/seurat_dimplot_colorby.R --seurat /tmp/pbmc.qs2 --colorby groups --reduction tsne --outdir /tmp/out
```

Plot scripts default to `--reduction umap`, so pass `--reduction tsne` when testing against the bundled objects.

## CLI script conventions

- `#!/usr/bin/env Rscript`, `pacman::p_load(optparse)`, `pdf(NULL)`, then a `# Purpose:` comment
- Options block between `##### Options #####` and `#####` banners; unpack `opts$x` into plain locals; validate before the (slow) object load; `print(opts)` under a `"\nArguments:"` message
- Load the analysis packages *after* argument parsing: `pacman::p_load(nvutils, sct2, Seurat, scCustomize, tidyverse)`
- Read objects with `ReadSeurat()` / write with `WriteSeurat()` — never `readRDS`/`qs2::qs_read` directly; both infer format from the extension so every script accepts `.rds` and `.qs2`
- `--seurat` is the input-object flag name
- Plot scripts write PNGs to `<outdir>/plots/` (`--outdir` default `.`); `ggsave(..., height = height, width = width, bg = "white")`
- End with `cat("\n\n")` then `devtools::session_info()`
- `chmod 755` new scripts, and add a section to the README's "Command-line scripts" list
