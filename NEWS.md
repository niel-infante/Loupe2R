# Loupe2R 0.5.0

## Breaking changes

- **An unrecognized `.cloupe` format version now aborts `cloupe_to_seurat()` with an error by default, instead of just warning.** Previously, a file reporting a format version outside the validated set emitted an R `warning()` and extraction proceeded anyway. `cloupe_to_seurat()` now stops before extracting anything unless called with the new `version_check = FALSE`, which restores the previous warn-and-proceed behavior (an R `warning()` is still raised) and puts the responsibility for verifying the result on you. This is driven by the same change in `cloupe_extract` 0.2.0 (raises `UnvalidatedFormatVersionError` in Python, which reticulate propagates as an R error) — upgrade `cloupe_extract` to pick it up.

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
