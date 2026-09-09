# Loupe2R - Claude Code Guide

## Project Overview

Loupe2R imports 10x Genomics `.cloupe` files (Visium HD, binned and cell-segmentation modes) into Seurat objects. All `.cloupe`-format parsing lives in the separate Python package [`loupe2py`](https://github.com/niel-infante/Loupe2Py) (this repo's sibling), called via `reticulate` — Loupe2R itself has no format-specific parsing code; it assembles a `Seurat` object from the generic files `loupe2py::extract_cloupe()` writes.

## Key Files

| File | Purpose |
|---|---|
| `R/cloupe_to_seurat.R` | Main entry point — calls `loupe2py` via reticulate, assembles the `Seurat` object |
| `R/combine_cloupe_bins.R` | Merges multiple single-resolution `cloupe_to_seurat()` outputs into one multi-assay object |
| `R/utils.R` | `.detect_mt_pattern()`, `.derive_sample_name()`, `.align_positions()` |
| `tests/testthat/test-integration-visium-hd.R` | Opt-in regression test against real paired `.cloupe`/SpaceRanger `outs/` data, gated by `LOUPE2R_TEST_DIR` |

## When a new SpaceRanger version is released

Because all `.cloupe`-format parsing lives in `loupe2py`, check **its** `CLAUDE.md` first — that's where the format-version guard, the `Matrices`-section synthetic-row handling, and the SpaceRanger-output-column checks belong.

Separately, re-run `test-integration-visium-hd.R` against a real sample from the new version: it reads official SpaceRanger output (`tissue_positions.parquet`, `scalefactors_json.json`, `filtered_feature_bc_matrix`) directly for cross-checking, so a new SpaceRanger version adding, renaming, or restructuring those columns could silently make this test's comparisons less thorough than they look, even if it still passes without error.

## Future work

**Build a real Seurat `FOV`/`Segmentation` object for cell-segmentation mode, not a centroid-only spatial assay.** `cloupe_to_seurat()` currently gives every cell-segmentation-mode cell a single `(x, y)` point (loupe2py 0.2.2 made that point exact, via the cell's true polygon area centroid — see `loupe2py`'s `CLAUDE.md`), then discards the actual cell shape. Seurat has a native object for the shape itself, confirmed directly against `SeuratObject`'s source (`github.com/satijalab/seurat-object`, `R/segmentation.R`), not guessed:

- **`CreateSegmentation(coords)`** — `coords` is a long-format data.frame with columns `cell`, `x`, `y`, **one row per polygon vertex**; it splits by `cell` and wraps each group into a polygon.
- **`CreateCentroids(coords)`** — same `cell`/`x`/`y` columns, but one row per cell (the point-based data we already have).
- **`CreateFOV(list(centroids = <Centroids obj>, segmentation = <Segmentation obj>))`** — combines both into one `FOV`, which is what `SpatialDimPlot(plot_segmentations = TRUE)` reads to draw the actual cell outlines instead of dots.

The needed input (`cell, x, y` per vertex) is a direct reshape of the `GeoJSON` polygon rings `loupe2py` already parses for the centroid fix — see that repo's `CLAUDE.md` for the extraction-side half of this (a new opt-in `cell_boundaries.csv` output). This side's work is: read that file, call the three functions above, attach the resulting `FOV` to the Seurat object in place of (or alongside) the current plain spatial assay for cell-segmentation-mode files specifically.

Not started. Real validation before considering this done: build a Seurat object from one of the real cell-segmented `.cloupe` files this project already has, and confirm `SpatialDimPlot(plot_segmentations = TRUE)` renders sane, non-degenerate polygons — not just that the object construction doesn't error.
