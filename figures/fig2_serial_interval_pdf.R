#!/usr/bin/env Rscript
#
# Fig 2 — Monthly forward serial-interval distribution (Claus Wilke ridgeline).
#   y axis in calendar order: Jan at the bottom -> Dec at the top.
#   (ggridges places the first factor level at the bottom, so
#    levels = month.abb gives Jan bottom, Dec top.)
#   alpha = 0.3, beta = 0.3 (stated in the manuscript caption, not on the figure).
#
#   Input : pipeline/processed/si_monthly_a0.3_b0.3.csv  (month, si_week, N, prob)
#
# Run: Rscript figures/fig2_serial_interval_pdf.R

suppressMessages({ library(readr); library(dplyr); library(ggplot2); library(ggridges) })

.args <- commandArgs(FALSE); .fa <- sub("^--file=", "", .args[grep("^--file=", .args)])
SCRIPT_DIR <- if (length(.fa)) dirname(normalizePath(.fa)) else getwd()
source(file.path(SCRIPT_DIR, "_common.R"))
OUT <- file.path(FIGDIR, "Fig2")

df <- read_csv(file.path(DATA, "si_monthly_a0.3_b0.3.csv"), show_col_types = FALSE) %>%
  mutate(mname = factor(month.abb[month], levels = month.abb))   # Jan bottom .. Dec top

p <- ggplot(df, aes(x = si_week, y = mname, height = prob, group = mname)) +
  geom_ridgeline(scale = 7, fill = SLATE, color = "white", linewidth = 0.3, min_height = 0) +
  scale_x_continuous("Serial interval (weeks)",
                     limits = c(0, 110), breaks = c(0, 30, 60, 90), expand = c(0, 0)) +
  scale_y_discrete(NULL, expand = expansion(mult = c(0.01, 0.12))) +
  coord_cartesian(clip = "off") +
  theme_ridges(font_size = 11, font_family = "Arial", grid = TRUE, center_axis_labels = TRUE) +
  theme(axis.text.y = element_text(vjust = 0), plot.margin = margin(6, 12, 6, 6))

save_fig(p, OUT, width = 5.4, height = 5.2)
