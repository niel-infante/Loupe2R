# Loupe2R - Claude Code Guide

## Project Overview

Loupe2R imports 10x Genomics `.cloupe` files (Visium HD, binned and cell-segmentation modes) into Seurat objects. All `.cloupe`-format parsing lives in [`cloupe_extract`](https://github.com/niel-infante/Loupe2Py/tree/main/cloupe_extract), a separate Python package that lives inside the `Loupe2Py` repo (this R repo's sibling) but is independently pip-installable — called via `reticulate`. `loupe2py` (the Python/AnnData tool) also depends on `cloupe_extract`, but Loupe2R does not depend on or import `loupe2py` itself, only `cloupe_extract` directly. Loupe2R has no format-specific parsing code of its own; it assembles a `Seurat` object from the generic files `cloupe_extract::extract_cloupe()` writes.

## Key Files

| File | Purpose |
|---|---|
| `R/cloupe_to_seurat.R` | Main entry point — calls `cloupe_extract` via reticulate, assembles the `Seurat` object |
| `R/combine_cloupe_bins.R` | Merges multiple single-resolution `cloupe_to_seurat()` outputs into one multi-assay object |
| `R/estimate_image_memory.R` | Calls `cloupe_extract`'s `estimate_image_memory()` via reticulate — scopes tissue-image memory cost across `px_per_bin` choices before committing to a real `cloupe_to_seurat()` call; see that package's own `CLAUDE.md` for the cheap-read mechanism |
| `R/utils.R` | `.detect_mt_pattern()`, `.derive_sample_name()`, `.align_positions()`, `.image_memory_table()` (reshapes `estimate_image_memory()`'s raw Python result into a data.frame) |
| `tests/testthat/test-integration-visium-hd.R` | Opt-in regression test against real paired `.cloupe`/SpaceRanger `outs/` data (binned mode), gated by `LOUPE2R_TEST_DIR` |
| `tests/testthat/test-integration-cellseg.R` | Opt-in regression test against real paired `.cloupe`/SpaceRanger cell-segmentation output, gated by `LOUPE2R_CELLSEG_TEST_DIR`. Uses `sf` for official polygon centroids (planar/GEOS, not spherical -- see the comment above `sf_use_s2(FALSE)` in the file). |

## When a new SpaceRanger version is released

**Validated so far: 3.0.0, 4.0.1, 4.1.0** (4.1.0 confirmed 2026-10-02 — see `cloupe_extract`'s `CLAUDE.md` for the full concordance numbers). Note that 4.1.0 emits the *same* internal `.cloupe` format versions as 4.0.1, so the format-version guard didn't fire and `_TESTED_VERSIONS` needed no change — meaning a clean `version_check` never tells you a SpaceRanger version was actually validated.

Because all `.cloupe`-format parsing lives in `cloupe_extract`, check **its** `CLAUDE.md` first — that's where the format-version guard, the `Matrices`-section synthetic-row handling, and the SpaceRanger-output-column checks belong. (`cloupe_extract` lives inside the `Loupe2Py` repo, under `cloupe_extract/`, alongside `loupe2py`'s own `CLAUDE.md` — check both if unsure which one covers something.)

Separately, re-run `test-integration-visium-hd.R` against a real sample from the new version: it reads official SpaceRanger output (`tissue_positions.parquet`, `scalefactors_json.json`, `filtered_feature_bc_matrix`) directly for cross-checking, so a new SpaceRanger version adding, renaming, or restructuring those columns could silently make this test's comparisons less thorough than they look, even if it still passes without error.

## Future work

**Build a real Seurat `FOV`/`Segmentation` object for cell-segmentation mode, not a centroid-only spatial assay.** `cloupe_to_seurat()` currently gives every cell-segmentation-mode cell a single `(x, y)` point (`cloupe_extract` 0.1.0 made that point exact, via the cell's true polygon area centroid — see `cloupe_extract`'s `CLAUDE.md`), then discards the actual cell shape. Seurat has a native object for the shape itself, confirmed directly against `SeuratObject`'s source (`github.com/satijalab/seurat-object`, `R/segmentation.R`), not guessed:

- **`CreateSegmentation(coords)`** — `coords` is a long-format data.frame with columns `cell`, `x`, `y`, **one row per polygon vertex**; it splits by `cell` and wraps each group into a polygon.
- **`CreateCentroids(coords)`** — same `cell`/`x`/`y` columns, but one row per cell (the point-based data we already have).
- **`CreateFOV(list(centroids = <Centroids obj>, segmentation = <Segmentation obj>))`** — combines both into one `FOV`, which is what `SpatialDimPlot(plot_segmentations = TRUE)` reads to draw the actual cell outlines instead of dots.

The needed input (`cell, x, y` per vertex) is a direct reshape of the `GeoJSON` polygon rings `cloupe_extract` already parses for the centroid fix — see that package's `CLAUDE.md` for the extraction-side half of this (a new opt-in `cell_boundaries.csv` output). This side's work is: read that file, call the three functions above, attach the resulting `FOV` to the Seurat object in place of (or alongside) the current plain spatial assay for cell-segmentation-mode files specifically.

Not started. Real validation before considering this done: build a Seurat object from one of the real cell-segmented `.cloupe` files this project already has, and confirm `SpatialDimPlot(plot_segmentations = TRUE)` renders sane, non-degenerate polygons — not just that the object construction doesn't error.

**Possible future thought: halve the R-side tissue-image memory footprint by reading the PNG as an integer array instead of double.** `cloupe_to_seurat()` loads the stitched image via `png::readPNG(img_path)`, whose default (`native = FALSE`) returns a `double`-typed array (8 bytes/channel/pixel) — this is what hits R's vector memory limit on large images (see `px_per_bin` above). Confirmed directly against the installed `Seurat` 5.5.1 namespace: `VisiumV1`/`VisiumV2`/`SliceImage`'s `image` slot is formally typed `"array"` with no validity function, and only checks that structural class — not `typeof`. Tested directly: a `nativeRaster` (from `readPNG(..., native = TRUE)`, ~1 packed int/pixel, by far the cheapest option) fails `is(x, "array")` and can't be used as-is; but a plain **integer**-typed array (same dims as today, values 0-255) passes `is(x, "array")` and successfully constructs a `VisiumV1` object — no code changes needed on Seurat's side. That would cut the double-array memory in half (4 bytes vs. 8 bytes per channel/pixel) on top of whatever `px_per_bin` already saves. Not verified: whether `SpatialDimPlot()`/`SpatialFeaturePlot()` render correctly off an integer-valued image array rather than the current double `[0,1]`-range one (R's `as.raster()` is documented to treat integer input as 0-255, so it should work, but this hasn't been test-rendered). Not scoped or started beyond this note.
