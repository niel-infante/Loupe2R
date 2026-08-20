# Loupe2R

Import 10x Genomics `.cloupe` files (including Visium HD) into [Seurat](https://satijalab.org/seurat/) objects.

`Loupe2R` extracts the count matrix, spatial coordinates, tissue image, UMAP embedding, and Space Ranger cluster labels from a `.cloupe` file and assembles a Seurat object from them, via a bundled Python extraction step run through [reticulate](https://rstudio.github.io/reticulate/). The actual `.cloupe` binary parsing is built directly on top of [`cellgeni/cloupe`](https://github.com/cellgeni/cloupe) — see [Credits](#credits).

**Read [Limitations & Risks](#limitations--risks) before using this on anything you plan to publish.** This package parses an undocumented, proprietary file format using an unofficial, reverse-engineered tool. It has been validated carefully on real data (see [Validation](#validation) below), but it is not a substitute for 10x Genomics' own SpaceRanger output when that's available.

## Installation

```r
# install.packages("remotes")
remotes::install_github("<your-org>/Loupe2R")
```

Loupe2R also requires, **separately**, a Python environment with `numpy`, `scipy`, and (optionally, for tissue images) `Pillow`, plus a local install of the [`cellgeni/cloupe`](https://github.com/cellgeni/cloupe) Python package (AGPL-3.0 — see [AGPL dependency](#agpl-dependency) below; **not** installed by Loupe2R):

```bash
pip install numpy scipy Pillow
pip install git+https://github.com/cellgeni/cloupe.git
```

By default Loupe2R looks for `cellgeni/cloupe` at `~/Tools/cloupe`. Point it elsewhere with `set_cloupe_path("/path/to/cloupe")` or `options(loupe2r.cloupe_pkg = "/path/to/cloupe")` in your `.Rprofile`.

## Quick start

```r
library(reticulate)
use_condaenv("base", required = TRUE)   # or any env with numpy/scipy/Pillow

library(Loupe2R)
srt <- cloupe_to_seurat("path/to/sample.cloupe")
Seurat::SpatialFeaturePlot(srt, features = "nCount_Spatial")
```

### Combining multiple Visium HD bin resolutions

A single `.cloupe` file only ever captures the one bin resolution that was open in Loupe Browser when it was saved. If you have separate `.cloupe` files for the same capture area at different resolutions, import each one and combine them:

```r
s008 <- cloupe_to_seurat("sample_008um.cloupe")
s016 <- cloupe_to_seurat("sample_016um.cloupe")
srt  <- combine_cloupe_bins(s008, s016)
Seurat::Assays(srt)  # "Spatial.008um" "Spatial.016um"
```

See `?combine_cloupe_bins` for the caveats this carries over from Seurat's own multi-resolution design (object-level cell identity tracks the first object passed in).

## What this does and does not do

- Extracts counts, spatial coordinates, the tissue image, UMAP/other embeddings, Space Ranger clusterings (prefixed `sr_`), and any user-created Loupe Browser cell tracks.
- Captures **one** bin resolution per `.cloupe` file — use `combine_cloupe_bins()` if you need more than one.
- Reflects whatever default, unfiltered pipeline output Loupe Browser was showing — not your own QC, filtering, or normalization choices.

## Limitations & Risks

**The `.cloupe` format is proprietary and undocumented.** It's a custom binary container (JSON header + byte-offset index, not SQLite) with no public spec. The parser this package depends on, [`cellgeni/cloupe`](https://github.com/cellgeni/cloupe), is an unofficial, reverse-engineered tool explicitly self-described upstream as a "PoC" (proof of concept), with no version-compatibility guarantees across Loupe Browser/CellRanger/SpaceRanger releases. A future 10x software update could change the format in a way this parser misreads. See [`check_format_version()`](#format-version-guard) below for how this package tries to catch that early.

**Loupe-derived annotations reflect default, unfiltered pipeline output, not your own analysis.** In [satijalab/seurat#9269](https://github.com/satijalab/seurat/issues/9269), a Seurat maintainer cautions that data exported from Loupe Browser reflects only the default parameters SpaceRanger/CellRanger ran with — no manual QC, filtering, or normalization decisions are captured. Treat `sr_*` clusterings and other Loupe-derived metadata as exploratory defaults, not a substitute for your own analysis.

**A `.cloupe` file cannot represent Visium HD's full multi-resolution structure.** Seurat's official path, `Seurat::Load10X_Spatial(bin.size = c(8, 16))`, loads multiple bin resolutions as separate assays from one SpaceRanger output directory. A `.cloupe` file only ever captures the single resolution that was open in Loupe Browser when it was saved; use `combine_cloupe_bins()` (above) if you have separate `.cloupe` files per resolution.

**Standard-Visium + cell-segmentation `.cloupe` files use a coarser spatial fallback.** Visium HD barcodes encode their own exact grid position (see [Validation](#validation)); non-HD cell-segmentation `.cloupe` files don't use that barcode convention, so their coordinates fall back to averaged cell-segment centroids, which is less precise.

**Prefer official SpaceRanger output when you have it.** If you have access to the SpaceRanger `outs/` directory for a sample, `Seurat::Load10X_Spatial()` on that is the authoritative path — use `Loupe2R` for `.cloupe`-only samples, or to cross-check.

### Format version guard

Loupe2R checks the `.cloupe` file's internal format-version fields against a known-tested set on every extraction. If a file reports a version outside that set, `cloupe_to_seurat()` emits an R `warning()` and still proceeds (an unrecognized version doesn't necessarily mean broken output, just untested) — but treat results from those files with extra scrutiny, and prefer cross-checking against official SpaceRanger output when possible. The detected versions are always stashed on the returned object via `Seurat::Misc(srt, "cloupe_format_info")`, regardless of whether they triggered a warning.

Validated against `.cloupe` files with container/index format `9.0.0`, run format `3.0.0`, matrix format `6.3.0`, and projection format `4.1.0` (checked across multiple independent Visium HD and standard-Visium samples). This is the maintenance model going forward: as files with other versions are validated, extend `_TESTED_VERSIONS` in `inst/python/cloupe_extract.py` — there's no way to preemptively support format versions that don't exist yet.

### AGPL dependency

Loupe2R itself is MIT-licensed. At runtime, it calls a separately-installed copy of [`cellgeni/cloupe`](https://github.com/cellgeni/cloupe) (AGPL-3.0) via reticulate/`sys.path`, at a path you configure. It is **not bundled, vendored, or redistributed** by this package — you install it yourself and point Loupe2R at it. This keeps Loupe2R's own MIT license clean of AGPL obligations; if you build anything that bundles or redistributes `cellgeni/cloupe` itself, that combined work is subject to AGPL-3.0.

## Credits

The hard part — reverse-engineering the proprietary `.cloupe` binary format at all (the header layout, byte-offset index block, matrix/projection encoding, and tiled image storage) — is not this package's work. It's [`cellgeni/cloupe`](https://github.com/cellgeni/cloupe), written by **Martin Prete** and **Nithin Mathew Joseph** of the Wellcome Sanger Institute's Cellular Genetics Informatics (cellgeni) team, AGPL-3.0 licensed. `Loupe2R` calls that parser directly and builds Visium HD-specific extraction on top of it — bin-size and tile-pyramid image handling, barcode-based grid-position parsing, and the format-version guard — plus the Seurat object assembly. None of that is possible without their groundwork; see [AGPL dependency](#agpl-dependency) above for how that dependency is used and licensed.

Also indebted to [Satija Lab's Seurat](https://satijalab.org/seurat/) and [10x Genomics' SpaceRanger](https://www.10xgenomics.com/support/software/space-ranger) for the reference outputs and object conventions this package's Visium HD support was validated against (see [Validation](#validation)).

## Validation

`cloupe_to_seurat()` was cross-validated against the official SpaceRanger output for the same Visium HD sample (8µm bins, 19,072 features × 472,752 barcodes). Results:

| Check | Result |
|---|---|
| Barcode overlap (cloupe vs. official filtered matrix) | 100% |
| `nCount_Spatial` correlation on shared barcodes | 1.000000 |
| `array_row`/`array_col` exact match vs. official `tissue_positions.parquet` | 100% (0 max difference) |
| `pxl_row_in_fullres`/`pxl_col_in_fullres` correlation | 1.000000 |
| `spot_diameter_fullres` vs. official `scalefactors_json.json` | exact match (1.721916) |

The `array_row`/`array_col` exact match is not a rounding coincidence: Visium HD barcodes encode their own grid position directly (e.g. `s_008um_00269_00526-1` → row 269, col 526), and `cloupe_extract.py` parses that directly rather than approximating it from pixel coordinates. This was in fact a real bug caught during development — a naive `round(pixel / bin_size_px)` reconstruction disagreed with the official grid by ~450 bins (a large, systematic offset, not boundary noise), because the stitched tissue image's pixel origin doesn't align with SpaceRanger's own coordinate reference frame. Parsing the barcode instead sidesteps that mismatch entirely. Non-HD barcodes (e.g. `cellid_...` in cell-segmentation `.cloupe` files) don't follow this convention and fall back to the pixel-based approximation instead.

Reproduce this validation yourself (see `tests/testthat/test-integration-visium-hd.R`) by setting `LOUPE2R_TEST_DIR` to a SpaceRanger Visium HD `outs/` directory containing `cloupe_008um.cloupe` and `binned_outputs/square_008um/`, and running the package's tests.

## License

MIT (this package) + file LICENSE. See [AGPL dependency](#agpl-dependency) above for the separately-installed `cellgeni/cloupe` requirement.
