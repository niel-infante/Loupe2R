# Loupe2R 0.7.0

## New features

- **`estimate_image_memory(cloupe_path, px_per_bin = list(NULL, 16, 8, 4, 2, 1))`**: scopes tissue-image memory cost across candidate `px_per_bin` resolutions before committing to a real `cloupe_to_seurat()` call. Calls `cloupe_extract`'s own `estimate_image_memory()` (0.4.0) via reticulate -- a cheap read that skips the count matrix entirely and never decodes tile pixel bytes, confirmed directly against the real CRC validation file at ~4.6s vs. `cloupe_to_seurat()`'s measured 667.7s. Returns a data.frame with `image_width`/`image_height`, `clamped_to_native`, and two memory estimates per resolution: `python_peak_gb` (the Pillow-side stitching peak) and `r_reload_peak_gb` (this package's own `png::readPNG()` reload peak -- the tighter constraint in practice, and the one that actually crashed on a real oversized image; see `cloupe_to_seurat()`'s `px_per_bin` documentation).

# Loupe2R 0.6.0

## Breaking changes

- **Tissue image reconstruction no longer defaults to native resolution.** `cloupe_to_seurat()` gained `px_per_bin` (default `4`, passed straight through to `cloupe_extract$extract_cloupe()`): the target tissue-image resolution, expressed as output pixels per finest-grid bin, rather than always reconstructing at full native resolution. On a real large Visium HD sample, unconditional native-resolution reconstruction was a confirmed crash: `png::readPNG()` expands the image to double precision and hit R's 30 GB vector memory limit outright. Pass `px_per_bin = NULL` to opt back into the old, only-ever behavior. `scale.factors$hires`/`$lowres` on the embedded `VisiumV2`/`VisiumV1` image are no longer always `1.0` -- they now reflect the real resize ratio; this package's coordinate handling already stored raw fullres coordinates and passed the scale factor through generically, so no other logic changed. See `cloupe_extract`'s own `NEWS.md` for the underlying fix.

## Testing

- New opt-in integration test, `test-integration-large-image.R` (gated by `LOUPE2R_LARGE_IMAGE_TEST_DIR`), is a direct regression test for the R-side crash above: runs `cloupe_to_seurat()` with the new default `px_per_bin` against the same large real sample that used to crash, and confirms the embedded image is genuinely smaller than native with a non-`1.0` scale factor.

# Loupe2R 0.5.0

## Breaking changes

- **An unrecognized `.cloupe` format version now aborts `cloupe_to_seurat()` with an error by default, instead of just warning.** Previously, a file reporting a format version outside the validated set emitted an R `warning()` and extraction proceeded anyway. `cloupe_to_seurat()` now stops before extracting anything unless called with the new `version_check = FALSE`, which restores the previous warn-and-proceed behavior (an R `warning()` is still raised) and puts the responsibility for verifying the result on you. This is driven by the same change in `cloupe_extract` 0.2.0 (raises `UnvalidatedFormatVersionError` in Python, which reticulate propagates as an R error) — upgrade `cloupe_extract` to pick it up.

## Testing

- New opt-in integration test, `test-integration-cellseg.R` (gated by `LOUPE2R_CELLSEG_TEST_DIR`), cross-validates `cloupe_to_seurat()` against real official SpaceRanger cell-segmentation output — barcode overlap, count concordance, and polygon-centroid position, the same methodology `test-integration-visium-hd.R` already applies to binned mode. Previously, cell-segmentation correctness for `Loupe2R` specifically was inferred from sharing `cloupe_extract` with `loupe2py` (which does have its own cell-segmentation test data), not independently confirmed. Run against the Human Kidney FFPE dataset used in the JBT paper draft: 100.0000% barcode overlap, 100.0000% count exact match, centroid distance 0.0000 px mean/median (148,056 cells) — matching the paper's reported numbers. Adds `sf` to `Suggests` (used for official polygon centroids; the test is skipped if not installed).

# Loupe2R 0.4.0

## Breaking changes

- `Loupe2R` now calls out to [`cloupe_extract`](https://github.com/niel-infante/Loupe2Py/tree/main/cloupe_extract) instead of `loupe2py`. The extraction core previously lived inside the `loupe2py` package itself — confusingly, since `loupe2py` also names the separate Python/AnnData sibling tool. It's now split out into its own, independently pip-installable package (still inside the `Loupe2Py` repo, under `cloupe_extract/`), which `Loupe2R` depends on directly and `loupe2py` depends on transitively.
- Install it with `pip install "cloupe_extract @ git+https://github.com/niel-infante/Loupe2Py.git#subdirectory=cloupe_extract"`. If you already have `loupe2py` installed from before this change, reinstall/upgrade it (`pip install --upgrade loupe2py`) to pick up the new `cloupe_extract` dependency, or install `cloupe_extract` directly — an old `loupe2py` install predating this split won't have it.

# Loupe2R 0.3.0

## License change

- `Loupe2R` is now licensed **AGPL-3.0-or-later** (previously MIT). It depends at runtime on `loupe2py`, which vendors AGPL-3.0-licensed code from [`cellgeni/cloupe`](https://github.com/cellgeni/cloupe); this release adopts that same license rather than drawing a technical line around the reticulate call. See the README's [AGPL dependency](https://github.com/niel-infante/Loupe2R#agpl-dependency) section for the reasoning.

# Loupe2R 0.2.0

## Breaking changes

- `Loupe2R` no longer bundles its own `.cloupe` extraction logic. It now depends on [`loupe2py`](https://github.com/niel-infante/Loupe2Py), a separate pip-installable Python package shared with the new squidpy/AnnData sibling tool. Install it with `pip install git+https://github.com/niel-infante/Loupe2Py.git` — no separate `cellgeni/cloupe` install or path configuration is needed anymore.
- `set_cloupe_path()` is removed, along with the `LOUPE2R_CLOUPE_PATH` env var and `options(loupe2r.cloupe_pkg = ...)`. There is no external path to configure now that `loupe2py` vendors its own pinned copy of the `.cloupe` parser.

## Other changes

- `cloupe_to_seurat()` now calls `reticulate::import("loupe2py")` directly instead of building a Python call via string interpolation, closing a latent string-injection risk from unescaped `.cloupe` file paths.

# Loupe2R 0.1.0

Initial release: `cloupe_to_seurat()` and `combine_cloupe_bins()`.
