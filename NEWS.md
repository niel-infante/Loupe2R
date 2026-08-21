# Loupe2R 0.2.0

## Breaking changes

- `Loupe2R` no longer bundles its own `.cloupe` extraction logic. It now depends on [`loupe2py`](https://github.com/niel-infante/Loupe2Py), a separate pip-installable Python package shared with the new squidpy/AnnData sibling tool. Install it with `pip install git+https://github.com/niel-infante/Loupe2Py.git` — no separate `cellgeni/cloupe` install or path configuration is needed anymore.
- `set_cloupe_path()` is removed, along with the `LOUPE2R_CLOUPE_PATH` env var and `options(loupe2r.cloupe_pkg = ...)`. There is no external path to configure now that `loupe2py` vendors its own pinned copy of the `.cloupe` parser.

## Other changes

- `cloupe_to_seurat()` now calls `reticulate::import("loupe2py")` directly instead of building a Python call via string interpolation, closing a latent string-injection risk from unescaped `.cloupe` file paths.

# Loupe2R 0.1.0

Initial release: `cloupe_to_seurat()` and `combine_cloupe_bins()`.
