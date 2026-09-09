# Loupe2R

Import 10x Genomics `.cloupe` files (including Visium HD) into [Seurat](https://satijalab.org/seurat/) objects.

Python/squidpy user? See the companion package, [**Loupe2Py**](https://github.com/niel-infante/Loupe2Py). Looking for the full function list (including internal, unexported functions)? See [`REFERENCE.md`](REFERENCE.md).

`Loupe2R` extracts the count matrix, spatial coordinates, tissue image, UMAP embedding, and Space Ranger cluster labels from a `.cloupe` file and assembles a Seurat object from them, via [reticulate](https://rstudio.github.io/reticulate/) calling out to [`loupe2py`](https://github.com/niel-infante/Loupe2Py) — a separate, pip-installable Python package that does the actual `.cloupe` binary parsing (shared with [Loupe2Py](https://github.com/niel-infante/Loupe2Py), the squidpy/AnnData sibling of this package). See [Credits](#credits).

**Read [Limitations & Risks](#limitations--risks) before using this on anything you plan to publish.** This package parses an undocumented, proprietary file format using an unofficial, reverse-engineered tool. It has been validated carefully on real data (see [Validation](#validation) below), but it is not a substitute for 10x Genomics' own SpaceRanger output when that's available.

## Installation

```r
# install.packages("remotes")
remotes::install_github("niel-infante/Loupe2R")
```

Loupe2R also requires, **separately**, a Python environment with [`loupe2py`](https://github.com/niel-infante/Loupe2Py) installed (**not** installed by Loupe2R itself):

```bash
pip install git+https://github.com/niel-infante/Loupe2Py.git
```

`loupe2py` pulls in `numpy`/`scipy`/`Pillow` itself and vendors its own copy of the `.cloupe` binary parser — no further Python-side setup, path configuration, or separate parser install is needed.

## Quick start

```r
library(reticulate)
use_condaenv("loupe2py", required = TRUE)   # any env with `pip install loupe2py`

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

**The `.cloupe` format is proprietary and undocumented.** It's a custom binary container (JSON header + byte-offset index, not SQLite) with no public spec. The parser this package ultimately depends on (vendored inside `loupe2py`; see [Credits](#credits)) is an unofficial, reverse-engineered tool with no version-compatibility guarantees across Loupe Browser/CellRanger/SpaceRanger releases. A future 10x software update could change the format in a way this parser misreads. See [Format version guard](#format-version-guard) below for how this package tries to catch that early.

**Loupe-derived annotations reflect default, unfiltered pipeline output, not your own analysis.** In [satijalab/seurat#9269](https://github.com/satijalab/seurat/issues/9269), a Seurat maintainer cautions that data exported from Loupe Browser reflects only the default parameters SpaceRanger/CellRanger ran with — no manual QC, filtering, or normalization decisions are captured. Treat `sr_*` clusterings and other Loupe-derived metadata as exploratory defaults, not a substitute for your own analysis.

**A `.cloupe` file cannot represent Visium HD's full multi-resolution structure.** Seurat's official path, `Seurat::Load10X_Spatial(bin.size = c(8, 16))`, loads multiple bin resolutions as separate assays from one SpaceRanger output directory. A `.cloupe` file only ever captures the single resolution that was open in Loupe Browser when it was saved; use `combine_cloupe_bins()` (above) if you have separate `.cloupe` files per resolution.

**Every `.cloupe` file's count matrix carried 2 bogus "genes" with inflated totals until `loupe2py` 0.1.1→0.2.1.** The `.cloupe` format itself appends synthetic `type_sum_<FeatureType>`/`genome_sum_<Reference>` rows to the `Matrices` section — each holding a barcode's *entire* total for that category, not a real gene's count. `loupe2py` wrote these out as if they were ordinary genes, which inflated `nCount_Spatial`/`percent.mt`-style metrics (3x, for the common single-feature-type/single-genome case) and added 2 meaningless entries to the gene list. Fixed in `loupe2py` 0.2.1 — see [its NEWS.md](https://github.com/niel-infante/Loupe2Py/blob/main/NEWS.md). Re-extract if you used an earlier version; the real per-gene counts themselves were always correct, only these 2 extra rows were the problem.

**Visium HD's cell-segmentation output uses the file's exact cell boundary data for position.** SpaceRanger 4.0+ can process a Visium HD run in two ways — binned (square bins at a chosen resolution) or cell-segmentation mode (transcripts assigned to individual cells via image-based nucleus/cell segmentation) — each producing its own separate `.cloupe` file; cell segmentation is a Visium HD analysis mode, not a separate non-HD product. Binned-mode barcodes encode their own exact grid position (see [Validation](#validation)); cell-segmentation-mode barcodes don't use that convention. Since `loupe2py` 0.2.2, position is computed as the true area centroid of each cell's embedded polygon boundary (`CellSegs.GeoJSON`) — confirmed coordinate-for-coordinate identical to SpaceRanger's own segmented output — which matches official ground truth exactly. Earlier versions averaged a coarser rectangular approximation of the cell mask instead (kept today only as a fallback for the rare barcode without polygon data): `loupe2py` 0.1.1 fixed a real bug in that averaging (`Centers` read with the wrong array layout, scrambling ~30-40% of cells onto a diagonal), and 0.2.2 replaced the averaging approach itself, which still had a small residual (~2 µm) even once correctly read; see [its NEWS.md](https://github.com/niel-infante/Loupe2Py/blob/main/NEWS.md). `array_row`/`array_col` are unconditionally `(0, 0)` for cell-segmentation-mode `.cloupe` files — there's no fixed grid for individual cells to round to — which doesn't affect plotting (that uses the real pixel coordinates) but is worth knowing if you read those columns directly.

**Prefer official SpaceRanger output when you have it.** If you have access to the SpaceRanger `outs/` directory for a sample, `Seurat::Load10X_Spatial()` on that is the authoritative path — use `Loupe2R` for `.cloupe`-only samples, or to cross-check.

### What a `.cloupe` file can't give you

Even for common cases, `.cloupe` doesn't carry everything in a full SpaceRanger output tree. Confirmed directly against real files:

| Not in `.cloupe` | What it's for |
|---|---|
| Raw (unfiltered) count matrix | Ambient-RNA correction, background-signal QC, custom re-filtering — `.cloupe`'s barcode set matches SpaceRanger's *filtered* matrix only |
| Molecule-level data (`molecule_info.h5`) | Per-UMI/per-read records for saturation analysis |
| Aligned reads (`possorted_genome_bam.bam`) | Variant calling, allele-specific expression, splice-junction/isoform analysis |
| Probe-level data (`probe_set.csv`, `raw_probe_bc_matrix.h5`) | Anything where gene-level aggregation is inappropriate (probe-based/FFPE chemistry) |
| `feature_slice.h5` | HD-specific rasterized per-gene expression image data |
| `barcode_mappings.parquet` | The exact HD bin-to-bin crosswalk (which 2 µm bins compose which 8 µm bin) — `combine_cloupe_bins()` merges resolutions as independent assays with no containment relationship between their bins, a practical but less precise substitute |
| Additional image products (`cytassist_image.tiff`, `aligned_fiducials.jpg`, `detected_tissue_image.jpg`) | Redoing tissue detection or checking fiducial-alignment quality |
| Full QC/metrics reporting (`web_summary.html`, `metrics_summary.csv`) | The `.cloupe` sections that would carry this (`Metrics`, `Analyses`) exist in the format but are confirmed empty in every real file checked — not partially present, genuinely absent |
| Differential expression / PCA tables | Confirmed absent, not just unextracted — the `Analyses` section contains only a null placeholder entry in every file checked |

### Format version guard

Loupe2R checks the `.cloupe` file's internal format-version fields against a known-tested set on every extraction. If a file reports a version outside that set, `cloupe_to_seurat()` emits an R `warning()` and still proceeds (an unrecognized version doesn't necessarily mean broken output, just untested) — but treat results from those files with extra scrutiny, and prefer cross-checking against official SpaceRanger output when possible. The detected versions are always stashed on the returned object via `Seurat::Misc(srt, "cloupe_format_info")`, regardless of whether they triggered a warning.

Validated against `.cloupe` files with container/index format `9.0.0`, run format `3.0.0`, matrix format `6.3.0`, and projection format `4.1.0` (checked across multiple independent Visium HD samples, in both binned and cell-segmentation modes). This is the maintenance model going forward: as files with other versions are validated, extend `_TESTED_VERSIONS` in [`loupe2py`'s `extract.py`](https://github.com/niel-infante/Loupe2Py/blob/main/src/loupe2py/extract.py) — there's no way to preemptively support format versions that don't exist yet.

### AGPL dependency

Loupe2R is licensed **AGPL-3.0-or-later** (see `LICENSE`). At runtime it calls a separately-installed Python package, `loupe2py`, via reticulate; `loupe2py` in turn vendors a pinned copy of [`cellgeni/cloupe`](https://github.com/cellgeni/cloupe)'s AGPL-3.0 `.cloupe` parser directly inside itself (see [Credits](#credits) and `loupe2py`'s own README). Loupe2R's functionality is inseparable from that AGPL-3.0 code in practice — there is no non-AGPL way to get a `.cloupe` file open — so Loupe2R adopts the same license rather than drawing a technical line around the reticulate call. This was a deliberate choice made when relicensing both packages together; MIT would likely have been legally defensible for Loupe2R specifically (it doesn't bundle or vendor the AGPL code itself), but AGPL-3.0 is the license that matches the spirit of the code both packages depend on.

## Credits

The hard part — reverse-engineering the proprietary `.cloupe` binary format at all (the header layout, byte-offset index block, matrix/projection encoding, and tiled image storage) — is not this package's work. It's [`cellgeni/cloupe`](https://github.com/cellgeni/cloupe), written by **Martin Prete** and **Nithin Mathew Joseph** of the Wellcome Sanger Institute's Cellular Genetics Informatics (cellgeni) team, AGPL-3.0 licensed, vendored inside [`loupe2py`](https://github.com/niel-infante/Loupe2Py) (which this package depends on). The Visium HD-specific extraction built on top of it — bin-size and tile-pyramid image handling, barcode-based grid-position parsing, and the format-version guard — also lives in `loupe2py`, shared with its squidpy/AnnData sibling; this package (`Loupe2R`) adds the Seurat object assembly on top. None of it is possible without cellgeni/cloupe's groundwork; see [AGPL dependency](#agpl-dependency) above for how that dependency is used and licensed.

Also indebted to [Satija Lab's Seurat](https://satijalab.org/seurat/) and [10x Genomics' SpaceRanger](https://www.10xgenomics.com/support/software/space-ranger) for the reference outputs and object conventions this package's Visium HD support was validated against (see [Validation](#validation)).

## Validation

`cloupe_to_seurat()` was cross-validated against the official SpaceRanger output for the same Visium HD sample (8µm bins, 19,072 features × 472,752 barcodes). Results:

| Check | Result |
|---|---|
| Barcode overlap (cloupe vs. official filtered matrix) | 100% |
| `nCount_Spatial` exact match vs. official (not just correlated) | 100% |
| `array_row`/`array_col` exact match vs. official `tissue_positions.parquet` | 100% (0 max difference) |
| `pxl_row_in_fullres`/`pxl_col_in_fullres` correlation | 1.000000 |
| `spot_diameter_fullres` vs. official `scalefactors_json.json` | exact match (1.721916) |

The `array_row`/`array_col` exact match is not a rounding coincidence: Visium HD barcodes encode their own grid position directly (e.g. `s_008um_00269_00526-1` → row 269, col 526), and `loupe2py`'s `parse_array_position()` parses that directly rather than approximating it from pixel coordinates. This was in fact a real bug caught during development — a naive `round(pixel / bin_size_px)` reconstruction disagreed with the official grid by ~450 bins (a large, systematic offset, not boundary noise), because the stitched tissue image's pixel origin doesn't align with SpaceRanger's own coordinate reference frame. Parsing the barcode instead sidesteps that mismatch entirely. Non-HD barcodes (e.g. `cellid_...` in cell-segmentation `.cloupe` files) don't follow this convention; for those, `array_row`/`array_col` are unconditionally `(0, 0)` rather than an approximation, since cell segments aren't laid out on any fixed grid to round to.

Reproduce this validation yourself (see `tests/testthat/test-integration-visium-hd.R`) by setting `LOUPE2R_TEST_DIR` to a SpaceRanger Visium HD `outs/` directory containing `cloupe_008um.cloupe` and `binned_outputs/square_008um/`, and running the package's tests.

## License

AGPL-3.0-or-later + file LICENSE. See [AGPL dependency](#agpl-dependency) above for why.
