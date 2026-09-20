# Loupe2R function reference

Every function in the package, public and internal. If you're calling something from outside the package, use only the **Public API** section — everything under **Internal** can change or disappear without notice and isn't part of any compatibility guarantee.

## Public API

Exported (`NAMESPACE`), documented (`?function_name` works), stable.

### `cloupe_to_seurat()`

```r
cloupe_to_seurat(
  cloupe_path,
  sample_name   = NULL,
  assay_name    = "Spatial",
  slice_name    = "slice1",
  include_image = TRUE,
  outdir        = NULL,
  keep_files    = FALSE,
  condaenv      = NULL,
  version_check = TRUE,
  px_per_bin    = 4
)
```

Imports one `.cloupe` file — Visium HD, in either binned or cell-segmentation mode — into a `Seurat` object. Calls out to [`cloupe_extract`](https://github.com/niel-infante/Loupe2Py/tree/main/cloupe_extract) via `reticulate` to do the actual `.cloupe` parsing, then assembles the count matrix, spatial coordinates, tissue image, UMAP/other embeddings, Space Ranger graph/k-means clusterings (stored as `sr_`-prefixed metadata columns), and any user-created Loupe Browser cell tracks into the returned object. Also stashes `bin_size_um` and `cloupe_format_info` via `Seurat::Misc()` for provenance, and computes `percent.mt`.

- `cloupe_path` — path to the `.cloupe` file.
- `sample_name` — stored in `orig.ident`; defaults to the `.cloupe` filename with its extension stripped.
- `assay_name`, `slice_name` — names for the Seurat assay and image slot.
- `include_image` — whether to reconstruct and embed the tissue image.
- `outdir`, `keep_files` — control where intermediate extracted files go; `NULL`/`FALSE` uses an auto-deleted temp directory.
- `condaenv` — reticulate conda environment to activate first, if not already configured.
- `version_check` — if `TRUE` (default), abort with an error when the file reports a `.cloupe` internal format version outside the validated set, before extracting anything. `FALSE` proceeds anyway (an R `warning()` is raised instead); you're then responsible for independently verifying the result.
- `px_per_bin` — target tissue-image resolution, as output pixels per finest-grid bin (or per nominal 10µm cell for cell-segmentation-mode files). Default `4` keeps the reconstructed image well clear of both a Python-side crash (Pillow's decompression-bomb guard) and an R-side one (`png::readPNG()` hitting R's vector memory limit) on large capture areas. Pass `NULL` for the old, only-ever behavior: full native resolution, no resize. Ignored when `include_image = FALSE`.

Returns a `Seurat` object.

### `combine_cloupe_bins()`

```r
combine_cloupe_bins(...)
```

Takes two or more `Seurat` objects, each returned by `cloupe_to_seurat()` on a different bin-resolution `.cloupe` file for the *same* tissue section, and merges them into one object with one assay and one image per resolution (`Spatial.008um`, `Spatial.016um`, ...), mirroring the structure `Seurat::Load10X_Spatial(bin.size = c(8, 16))` produces from a full SpaceRanger output directory. Errors if two inputs resolve to the same bin size; warns if bin counts across resolutions look inconsistent for the same tissue area (a guard against accidentally combining different samples — nothing in a single `.cloupe` file identifies which capture area it came from). The combined object's own cell identity (`Cells()`/`ncol()`, `meta.data`, idents) stays pinned to whichever object is passed first; index into a specific assay for resolution-specific work.

Returns a `Seurat` object.

## Internal

Not exported, not documented with a `?` help page, no compatibility guarantee — listed here for completeness and because reading the source is easier when you know what's already been built.

| Function | File | What it does |
|---|---|---|
| `.detect_mt_pattern(gene_names)` | `R/utils.R` | Returns `"^MT-"` if any gene name starts with `MT-` (human convention), else `"^mt-"` (mouse convention) — even if nothing actually matches, in which case `percent.mt` ends up 0. Used by `cloupe_to_seurat()`. |
| `.derive_sample_name(cloupe_path)` | `R/utils.R` | Strips the directory and final extension from a `.cloupe` path to derive the default `sample_name`. |
| `.align_positions(pos, cell_names)` | `R/utils.R` | Reindexes a `data.frame` with a `barcode` column to a given cell-name order — used to align `tissue_positions.csv` rows to the count matrix's column order. |
| `.warn_if_bin_counts_inconsistent(objs, bin_um)` | `R/combine_cloupe_bins.R` | Called by `combine_cloupe_bins()`; warns if the bin-count ratio between two resolutions doesn't match what's expected for the same tissue area. |
| `.onAttach(libname, pkgname)` | `R/zzz.R` | R package lifecycle hook, run automatically once when the package is loaded via `library(Loupe2R)` — prints the startup caution about the reverse-engineered `.cloupe` format. Not something you call directly. |
