# Install R packages needed for local manuscript figure-making.


if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

pacman::p_load(
  tidyverse,
  cowplot,
  effsize,
  ggdist,
  ggpubr,
  ggseg,
  ggsegSchaefer,
  gridExtra,
  gratia,
  knitr,
  kableExtra,
  mgcv,
  paletteer,
  pammtools,
  purrr,
  RColorBrewer,
  rjson,
  scales,
  PNWColors,
  pals,
  ggthemes
)

#################
# Pinned plotting packages for the cortical surface figures
#################
# The figure code depends on ggsegSchaefer v2.0.4 (2022 release). Newer
# ggsegSchaefer releases are a rewrite that is incompatible with this code,
# so install the versions below rather than the latest from CRAN/GitHub.
# Note: ggseg requires the sf package, which needs GDAL, GEOS, and PROJ
# installed at the system level (e.g., via Homebrew on macOS).

# macOS only: if compiling packages from source fails with missing header
# errors, point R to the Command Line Tools SDK
if (Sys.info()[["sysname"]] == "Darwin") {
  Sys.setenv(SDKROOT = system("xcrun --show-sdk-path", intern = TRUE))
}

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

# pin the plotting stack (install in this order so ggplot2 stays at 3.5.2)
remotes::install_version("ggplot2", version = "3.5.2", upgrade = "never")
remotes::install_version("ggseg",   version = "1.6.5", upgrade = "never")
remotes::install_version("ggseg3d", version = "1.6.3", upgrade = "never")

# install ggsegSchaefer v2.0.4 from GitHub
remotes::install_github(
  "ggsegverse/ggsegSchaefer@v2.0.4",
  dependencies = FALSE,
  upgrade = "never",
  build_vignettes = FALSE
)

# restart R, then check versions: should print 2.0.4, 1.6.5, 1.6.3, 3.5.2
sapply(c("ggsegSchaefer", "ggseg", "ggseg3d", "ggplot2"),
       function(p) as.character(packageVersion(p)))