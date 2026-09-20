#' Import a 10x Genomics .cloupe file into a Seurat object
#'
#' Calls the Python package \href{https://github.com/niel-infante/Loupe2Py/tree/main/cloupe_extract}{cloupe_extract}
#' (the shared .cloupe extraction core also used by loupe2py, this package's
#' Python/AnnData sibling) via reticulate to read the .cloupe binary, then
#' assembles a Seurat object with the count matrix, spatial coordinates,
#' tissue image, UMAP embedding, and Spaceranger cluster labels.
#'
#' @param cloupe_path   Path to the .cloupe file.
#' @param sample_name   Sample name stored in \code{orig.ident}. NULL (default)
#'                      derives the name from the .cloupe filename.
#' @param assay_name    Name for the Seurat assay (default "Spatial").
#' @param slice_name    Name for the image slot (default "slice1").
#' @param include_image Whether to extract and embed the tissue image.
#' @param outdir        Directory for intermediate extracted files.
#'                      NULL = auto temp dir (deleted on return).
#' @param keep_files    If TRUE and outdir is NULL, keep intermediate files.
#'                      Useful for debugging.
#' @param condaenv      Name of conda environment to activate before extraction
#'                      (NULL = use current reticulate Python). Set this to
#'                      whichever environment has \code{cloupe_extract} installed
#'                      if you have not otherwise configured reticulate.
#' @param version_check If TRUE (default), abort with an error when the file's
#'                      internal .cloupe format version(s) haven't been
#'                      validated against real paired SpaceRanger output.
#'                      Set to FALSE to proceed anyway (a warning is still
#'                      raised and the detected versions are still stashed via
#'                      \code{Seurat::Misc(srt, "cloupe_format_info")}) -- you
#'                      are then responsible for independently verifying the
#'                      result before trusting it.
#' @param px_per_bin    Target tissue-image resolution, as output pixels per
#'                      finest-grid bin (or per nominal 10um cell for
#'                      cell-segmentation-mode files). Default 4 keeps the
#'                      reconstructed image well clear of both a Python-side
#'                      crash (Pillow's decompression-bomb guard tripping on
#'                      reload) and an R-side one (\code{png::readPNG()}
#'                      expanding a very large image to double precision and
#'                      hitting R's vector memory limit) -- full
#'                      native-resolution reconstruction can be several
#'                      billion pixels for a large capture area and is a real
#'                      crash risk on both sides. Pass NULL for the old,
#'                      only-ever behavior: full native resolution, no
#'                      resize. Ignored when \code{include_image = FALSE}.
#'
#' @return A Seurat object with metadata columns \code{orig.ident} and
#'   \code{percent.mt} in addition to the standard \code{nFeature_Spatial}
#'   and \code{nCount_Spatial} fields.
#'
#' @examples
#' \dontrun{
#' library(reticulate)
#' use_condaenv("loupe2py", required = TRUE)  # any env with `cloupe_extract` installed
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
  condaenv      = NULL,
  version_check = TRUE,
  px_per_bin    = 4
) {
  if (!is.null(condaenv))
    reticulate::use_condaenv(condaenv, required = TRUE)

  cloupe_path <- normalizePath(path.expand(cloupe_path), mustWork = TRUE)

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
  cloupe_extract <- tryCatch(
    reticulate::import("cloupe_extract"),
    error = function(e) stop(
      "Failed to import cloupe_extract. Install it with:\n",
      "  pip install \"cloupe_extract @ git+https://github.com/niel-infante/Loupe2Py.git#subdirectory=cloupe_extract\"\n",
      "into the Python environment reticulate will use (numpy, scipy, and\n",
      "Pillow are pulled in automatically as cloupe_extract dependencies; no\n",
      "separate cloupe-parser install or path configuration is needed). This\n",
      "is also installed automatically if you instead `pip install loupe2py`\n",
      "(the Python/AnnData sibling package, which depends on cloupe_extract).\n\n",
      "Original error: ", conditionMessage(e)
    )
  )

  message(sprintf("Extracting from %s (may take several minutes)...",
                  basename(cloupe_path)))
  # version_check = TRUE (default) raises a Python UnvalidatedFormatVersionError
  # here -- which reticulate propagates as an R error, aborting this function
  # immediately -- if the file reports a .cloupe format version outside the
  # validated set. Pass version_check = FALSE to proceed anyway; see the
  # warning block below for what happens in that case.
  # px_per_bin = NULL is passed through as-is; reticulate converts R NULL to
  # Python None automatically, which extract_cloupe() treats as an explicit
  # opt-out into native-resolution reconstruction (see its own docstring).
  cloupe_extract$extract_cloupe(
    cloupe_path, outdir, include_image = include_image, version_check = version_check,
    px_per_bin = px_per_bin
  )

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
  # Format-version provenance. With version_check = TRUE (default), the
  # extract_cloupe() call above already aborted this function entirely if
  # anything was unvalidated -- so warnings here only appear when the caller
  # explicitly passed version_check = FALSE and accepted that risk. Either
  # way, stash the detected versions + bin size on the object for
  # methods-section provenance.
  # ------------------------------------------------------------------
  fmt_info_path <- file.path(outdir, "format_info.json")
  if (file.exists(fmt_info_path)) {
    fmt_info <- jsonlite::read_json(fmt_info_path, simplifyVector = TRUE)
    if (length(fmt_info$warnings) > 0) {
      warning(
        "Unvalidated .cloupe format version(s) detected (proceeding because ",
        "version_check = FALSE):\n  ",
        paste(fmt_info$warnings, collapse = "\n  "),
        "\nExtraction proceeded, but you are responsible for independently ",
        "verifying results with extra scrutiny. See the package README for ",
        "the list of validated format versions.",
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
    # spot = spot_d * hires_sf is the general, correct formula for that:
    # spot_d is always in fullres pixel units, and hires_sf (now genuinely
    # < 1.0 whenever px_per_bin downsamples the image, see extract.py) rescales
    # it down to the embedded raster's own pixel units for plotting.
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
