#!/usr/bin/env Rscript
#
# S1 Fig — Distribution of I_M + T_MH (mosquito latent period + M->H transmission
#   delay, in days) across the (alpha, beta) sensitivity grid (3x3).
#   9 combinations (alpha, beta) in {0.2, 0.3, 0.4}^2. Each panel: day 0-60
#   probability mass (normalised PDF, sums to 1) as a filled area (no annotations).
#   h_HM = 1, t_1 in 2013-2025. Slate fill (#5B7DA3), matching Fig 2.
#   Grid layout matches S2/S3 Fig (facet_grid beta-rows x alpha-cols).
#
#   Input : pipeline/processed/im_tmh_hist.csv  (alpha, beta, day, prob)
#
# Run: Rscript figures/figS1_IM_plus_TMH_pdf.R

suppressMessages({ library(readr); library(dplyr); library(ggplot2) })

.args <- commandArgs(FALSE); .fa <- sub("^--file=", "", .args[grep("^--file=", .args)])
SCRIPT_DIR <- if (length(.fa)) dirname(normalizePath(.fa)) else getwd()
source(file.path(SCRIPT_DIR, "_common.R"))
OUT <- file.path(FIGDIR, "S1_Fig"); FILL <- SLATE

h <- read_csv(file.path(DATA, "im_tmh_hist.csv"), show_col_types = FALSE) %>%
  mutate(alpha = factor(alpha), beta = factor(beta))

XMAX <- 60
XLAB <- expression(italic(I)[M] + italic(T)[MH] ~ "(days)")

p <- ggplot(h, aes(day, prob)) +
  geom_area(fill = FILL, colour = NA) +
  facet_grid(beta ~ alpha,
             labeller = label_bquote(cols = alpha == .(as.character(alpha)),
                                     rows = beta == .(as.character(beta)))) +
  scale_x_continuous(XLAB, breaks = c(0, 20, 40, 60), limits = c(0, XMAX), expand = c(0.01, 0)) +
  scale_y_continuous("Probability", breaks = c(0, 0.05, 0.10),
                     expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(ylim = c(0, 0.13)) +
  theme_minimal_hgrid(font_size = 10, font_family = "Arial") +
  theme(panel.spacing = unit(8, "pt"), axis.text = element_text(size = 8),
        axis.line.x = element_line(linewidth = 0.4),
        axis.title.x = element_text(margin = margin(t = 8)),
        axis.title.y = element_text(margin = margin(r = 8)),
        aspect.ratio = 3.5 / 5.8)

save_fig(p, OUT, width = 7.4, height = 5.2)
