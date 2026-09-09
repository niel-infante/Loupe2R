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
