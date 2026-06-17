#
# _common.R — shared setup for all figure scripts (Claus Wilke style, PLOS ONE spec)
#
#   - Palette : Okabe-Ito (Wilke's default colour-blind-safe palette)
#   - Font    : Arial (PLOS ONE requirement). With XQuartz installed, the cairo
#               devices (cairo_ps / cairo_pdf) subset-embed the real Arial font,
#               so text stays as selectable/searchable vectors.
#   - Paths   : the figure inputs are the CSVs produced by the Julia pipeline in
#               pipeline/processed/ (run pipeline/run_pipeline.sh first); the raw
#               weekly case series is read from data/raw_inputs/. Outputs go to
#               output/fig/. Paths are resolved relative to this file so the
#               scripts run from any working directory.
#   - Export  : save_fig() writes a PLOS-compliant EPS (vector, Arial embedded)
#               and TIFF (RGB 8-bit, no alpha, LZW, 300 dpi) by default.
#
# This file is sourced by every figures/figNN_*.R script.

suppressMessages({ library(ggplot2); library(cowplot) })

# ── path resolution ──────────────────────────────────────────────────────────
.args <- commandArgs(FALSE)
.fa   <- sub("^--file=", "", .args[grep("^--file=", .args)])
SCRIPT_DIR <- if (length(.fa)) dirname(normalizePath(.fa)) else getwd()
ROOT   <- normalizePath(file.path(SCRIPT_DIR, ".."))
DATA   <- file.path(ROOT, "pipeline", "processed")   # figure inputs (produced by the pipeline)
RAW    <- file.path(ROOT, "data", "raw_inputs")
FIGDIR <- file.path(ROOT, "output", "fig")

# ── colours (Okabe-Ito) ──────────────────────────────────────────────────────
OKABE <- c(orange = "#E69F00", skyblue = "#56B4E9", green = "#009E73",
           yellow = "#F0E442", blue = "#0072B2", vermillion = "#D55E00",
           purple = "#CC79A7", black = "#000000")

# Recurring project colours (kept consistent across figures):
#   serial-interval / probability distributions = slate  #5B7DA3
#   same-year / short IIP                        = steel  #3D6A93
#   cross-year / long IIP                        = terracotta #C56B4A
SLATE      <- "#5B7DA3"
STEEL      <- "#3D6A93"
TERRACOTTA <- "#C56B4A"

# ── export wrapper (vector EPS + raster TIFF, physical size preserved) ────────
# Mirrors cowplot::save_plot's spirit: specify physical size (inches) directly
# for reproducibility.
#
#   EPS  : cairo_ps with Arial embedded (PLOS-preferred vector text).
#   TIFF : ragg::agg_tiff opened directly (NOT via ggsave). ggsave-to-TIFF
#          produces RGBA (4 channels) → an alpha channel that violates PLOS;
#          opening the device directly with a white background flattens to RGB.
save_fig <- function(plot, stem, width, height,
                     eps = TRUE, tiff = TRUE, pdf = FALSE, png = FALSE) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  if (eps)  ggplot2::ggsave(paste0(stem, ".eps"), plot, device = cairo_ps,
                            width = width, height = height,
                            family = "Arial", fallback_resolution = 1200)
  if (tiff) {
    ragg::agg_tiff(paste0(stem, ".tiff"), width = width, height = height,
                   units = "in", res = 300, compression = "lzw",
                   background = "white")
    print(plot); grDevices::dev.off()
  }
  if (pdf)  ggplot2::ggsave(paste0(stem, ".pdf"), plot, device = cairo_pdf,
                            width = width, height = height, family = "Arial")
  if (png)  ggplot2::ggsave(paste0(stem, ".png"), plot, device = ragg::agg_png,
                            width = width, height = height, dpi = 400,
                            background = "white")
  cat("Saved:", basename(stem),
      if (eps) ".eps" else "", if (tiff) ".tiff" else "",
      if (pdf) ".pdf" else "", if (png) ".png" else "", "\n")
}
