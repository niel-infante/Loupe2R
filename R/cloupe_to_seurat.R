#' Import a 10x Genomics .cloupe file into a Seurat object
#'
#' Calls a Python extraction step via reticulate to read the .cloupe binary,
#' then assembles a Seurat object with the count matrix, spatial coordinates,
#' tissue image, UMAP embedding, and Spaceranger cluster labels.
#'
#' @param cloupe_path   Path to the .cloupe file.
#' @param sample_name   Sample name stored in \code{orig.ident}. NULL (default)
#'                      derives the name from the .cloupe filename.
#' @param assay_name    Name for the Seurat assay (default "Spatial").
#' @param slice_name    Name for the image slot (default "slice1").
#' @param include_image Whether to extract and embed the tissue image.
#'                      Requires Pillow in the active Python environment.
#' @param outdir        Directory for intermediate extracted files.
#'                      NULL = auto temp dir (deleted on return).
#' @param keep_files    If TRUE and outdir is NULL, keep intermediate files.
#'                      Useful for debugging.
#' @param condaenv      Name of conda environment to activate before extraction
#'                      (NULL = use current reticulate Python). Set to "base"
#'                      if you have not configured reticulate yet.
#'
#' @return A Seurat object with metadata columns \code{orig.ident} and
#'   \code{percent.mt} in addition to the standard \code{nFeature_Spatial}
#'   and \code{nCount_Spatial} fields.
#'
#' @examples
#' \dontrun{
#' library(reticulate)
#' use_condaenv("base", required = TRUE)   # base conda has numpy/scipy/Pillow
#'
#' srt <- cloupe_to_seurat("path/to/sample.cloupe")
#' SpatialFeaturePlot(srt, features = "nCount_Spatial")
#' }
#'
#' @export
#' @importFrom methods new
#' @importFrom utils read.csv
cloupe_to_seurat <- function(
  cloupe_path,
  sample_name   = NULL,
  assay_name    = "Spatial",
  slice_name    = "slice1",
  include_image = TRUE,
  outdir        = NULL,
  keep_files    = FALSE,
  condaenv      = NULL
) {
  if (!is.null(condaenv))
    reticulate::use_condaenv(condaenv, required = TRUE)

  cloupe_path <- normalizePath(path.expand(cloupe_path), mustWork = TRUE)

  # Locate the Python extraction helper bundled with the package
  py_helper <- system.file("python", "cloupe_extract.py", package = "Loupe2R")
  if (!nzchar(py_helper))
    stop("cloupe_extract.py not found in Loupe2R installation. Reinstall the package.")

  # Set up output directory
  cleanup <- is.null(outdir) && !keep_files
  if (is.null(outdir)) {
    outdir <- tempfile(pattern = "cloupe_")
    dir.create(outdir)
  }
  if (cleanup) on.exit(unlink(outdir, recursive = TRUE), add = TRUE)

  # ------------------------------------------------------------------
  # Python extraction
  # ------------------------------------------------------------------
  # Fail fast with an actionable message if the cloupe package isn't where
  # it's expected, rather than letting reticulate surface a cryptic
  # ImportError from deep inside Python.
  cloupe_pkg <- .resolve_cloupe_pkg()
  Sys.setenv(LOUPE2R_CLOUPE_PATH = cloupe_pkg)

  # Add the directory containing cloupe_extract.py to Python's path, then import
  helper_dir <- dirname(py_helper)
  reticulate::py_run_string(sprintf(
    'import sys as _sys\n_h = "%s"\nif _h not in _sys.path: _sys.path.insert(0, _h)',
    helper_dir
  ))

  tryCatch(
    reticulate::py_run_string("import cloupe_extract as _ce"),
    error = function(e) stop(
      "Failed to import cloupe_extract. Check that the Python environment has\n",
      "numpy, scipy, Pillow installed, and that the cloupe package is at:\n  ",
      cloupe_pkg, "\n\nOriginal error: ", conditionMessage(e)
    )
  )

  message(sprintf("Extracting from %s (may take several minutes)...",
                  basename(cloupe_path)))
  reticulate::py_run_string(sprintf(
    '_ce.extract_cloupe("%s", "%s", include_image=%s)',
    cloupe_path,
    outdir,
    if (include_image) "True" else "False"
  ))

  # ------------------------------------------------------------------
  # Build Seurat object from extracted files
  # ------------------------------------------------------------------
  message("Building Seurat object...")

  mat <- Seurat::ReadMtx(
    mtx            = file.path(outdir, "matrix.mtx.gz"),
    cells          = file.path(outdir, "barcodes.tsv.gz"),
    features       = file.path(outdir, "features.tsv.gz"),
    feature.column = 2,
    cell.column    = 1,
    skip.cell      = 0,
    skip.feature   = 0
  )

  sf  <- jsonlite::read_json(file.path(outdir, "scalefactors_json.json"))
  pos <- read.csv(file.path(outdir, "tissue_positions.csv"),
                  stringsAsFactors = FALSE)
  pos <- .align_positions(pos, colnames(mat))

  meta <- data.frame(
    in_tissue          = as.integer(pos$in_tissue),
    array_row          = pos$array_row,
    array_col          = pos$array_col,
    pxl_row_in_fullres = pos$pxl_row_in_fullres,
    pxl_col_in_fullres = pos$pxl_col_in_fullres,
    row.names          = colnames(mat)
  )

  srt <- Seurat::CreateSeuratObject(counts = mat, assay = assay_name, meta.data = meta)

  sname <- if (!is.null(sample_name)) sample_name else .derive_sample_name(cloupe_path)
  srt$orig.ident <- sname
  Seurat::Idents(srt) <- "orig.ident"

  mt_pat <- .detect_mt_pattern(rownames(srt))
  srt[["percent.mt"]] <- Seurat::PercentageFeatureSet(srt, pattern = mt_pat,
                                                       assay = assay_name)
  srt$percent.mt[is.na(srt$percent.mt)] <- 0

  # ------------------------------------------------------------------
  # Format-version provenance: surface as a warning if anything in the file
  # is outside the range of .cloupe versions this package has been
  # validated against, and stash the detected versions + bin size on the
  # object either way (useful for debugging / methods-section provenance).
  # ------------------------------------------------------------------
  fmt_info_path <- file.path(outdir, "format_info.json")
  if (file.exists(fmt_info_path)) {
    fmt_info <- jsonlite::read_json(fmt_info_path, simplifyVector = TRUE)
    if (length(fmt_info$warnings) > 0) {
      warning(
        "Unvalidated .cloupe format version(s) detected:\n  ",
        paste(fmt_info$warnings, collapse = "\n  "),
        "\nExtraction proceeded, but treat results with extra scrutiny. ",
        "See the package README for the list of validated format versions.",
        call. = FALSE
      )
    }
    Seurat::Misc(srt, "cloupe_format_info") <- fmt_info$versions
  }
  if (!is.null(sf[["bin_size_um"]])) {
    Seurat::Misc(srt, "bin_size_um") <- sf[["bin_size_um"]]
  }

  # ------------------------------------------------------------------
  # Embed tissue image
  # ------------------------------------------------------------------
  img_path <- file.path(outdir, "tissue_hires_image.png")

  if (include_image && file.exists(img_path)) {
    message("Loading tissue image...")
    img_arr  <- png::readPNG(img_path)
    hires_sf <- sf[["tissue_hires_scalef"]]
    spot_d   <- sf[["spot_diameter_fullres"]]

    # Store raw, unscaled fullres pixel coordinates (Seurat's own convention
    # per GetTissueCoordinates()/Read10X_Coordinates() -- scaling to whatever
    # raster is displayed happens at plot time, not at construction time).
    # hires_sf is always 1.0 here today (see the comment in cloupe_extract.py
    # on why), so this was previously a harmless no-op, but multiplying by it
    # here would double-scale coordinates if that ever changed.
    coords_df <- data.frame(
      tissue   = 1L,
      row      = pos$array_row,
      col      = pos$array_col,
      imagerow = pos$pxl_row_in_fullres,
      imagecol = pos$pxl_col_in_fullres,
      row.names = rownames(pos)
    )
    scale_facs <- Seurat::scalefactors(
      spot     = spot_d * hires_sf,
      fiducial = spot_d * hires_sf,
      hires    = hires_sf,
      lowres   = sf[["tissue_lowres_scalef"]]
    )

    image_obj <- tryCatch(
      new("VisiumV2", image = img_arr, scale.factors = scale_facs,
          coordinates = coords_df, assay = assay_name,
          key = paste0(slice_name, "_")),
      error = function(e)
        new("VisiumV1", image = img_arr, scale.factors = scale_facs,
            coordinates = coords_df, assay = assay_name,
            key = paste0(slice_name, "_"))
    )

    # Modern Seurat/SeuratObject (>= 5.4) has no explicit "default image"
    # API to set -- SpatialFeaturePlot() and friends resolve which image(s)
    # to use per-assay automatically (images = NULL -> Images(object,
    # assay = DefaultAssay(object))), so no further action is needed here.
    srt[[slice_name]] <- image_obj

  } else if (include_image) {
    warning("No tissue image found. Install Pillow: pip install Pillow")
  }

  # ------------------------------------------------------------------
  # Add non-spatial projections (UMAP, tSNE, ...) as DimReduc objects
  # ------------------------------------------------------------------
  proj_path <- file.path(outdir, "projections.csv")
  if (file.exists(proj_path)) {
    proj <- read.csv(proj_path, row.names = 1, stringsAsFactors = FALSE)
    proj <- proj[colnames(srt), , drop = FALSE]
    proj_names <- unique(sub("_[0-9]+$", "", colnames(proj)))
    for (pn in proj_names) {
      pn_cols <- grep(paste0("^", pn, "_[0-9]+$"), colnames(proj), value = TRUE)
      embed   <- as.matrix(proj[, pn_cols, drop = FALSE])
      key     <- paste0(tolower(pn), "_")
      colnames(embed) <- paste0(key, seq_len(ncol(embed)))
      srt[[tolower(pn)]] <- Seurat::CreateDimReducObject(
        embeddings = embed, key = key, assay = assay_name
      )
    }
    message(sprintf("Added reductions: %s", paste(tolower(proj_names), collapse = ", ")))
  }

  # ------------------------------------------------------------------
  # Add Spaceranger cluster labels
  # ------------------------------------------------------------------
  cl_path <- file.path(outdir, "clusterings.csv")
  if (file.exists(cl_path)) {
    cl_data <- read.csv(cl_path, row.names = 1, stringsAsFactors = FALSE)
    cl_data <- cl_data[colnames(srt), , drop = FALSE]
    for (cn in colnames(cl_data))
      srt[[paste0("sr_", make.names(cn))]] <- as.factor(cl_data[[cn]])
    message(sprintf("Added %d Spaceranger clustering(s) (prefix 'sr_')", ncol(cl_data)))
  }

  # ------------------------------------------------------------------
  # Add user cell tracks
  # ------------------------------------------------------------------
  ct_path <- file.path(outdir, "celltracks.csv")
  if (file.exists(ct_path)) {
    ct_data <- read.csv(ct_path, row.names = 1, stringsAsFactors = FALSE)
    ct_data <- ct_data[colnames(srt), , drop = FALSE]
    for (cn in colnames(ct_data))
      srt[[make.names(cn)]] <- ct_data[[cn]]
    message(sprintf("Added %d user cell track(s)", ncol(ct_data)))
  }

  message(sprintf(
    "Done!  %d genes x %d spots  |  assay: %s  |  image: %s",
    nrow(srt), ncol(srt), assay_name,
    if (slice_name %in% Seurat::Images(srt)) slice_name else "none"
  ))
  return(srt)
}
