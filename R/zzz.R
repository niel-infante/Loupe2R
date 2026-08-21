.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "Loupe2R imports .cloupe files via loupe2py, which wraps an unofficial,\n",
    "reverse-engineered parser for .cloupe's undocumented binary format. A\n",
    ".cloupe file captures only ONE bin resolution (not the multi-resolution\n",
    "structure Seurat::Load10X_Spatial(bin.size=c(8,16)) gives you) -- use\n",
    "combine_cloupe_bins() if you have separate files per resolution. See the\n",
    "package README/vignette before trusting results beyond exploration.\n",
    "(Silence with suppressPackageStartupMessages().)"
  )
}
