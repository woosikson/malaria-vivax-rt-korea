#!/usr/bin/env Rscript
#
# S4 Fig — Monthly mean short-IIP / long-IIP cases (alpha=beta=0.3), two solid lines.
#   short_iip = N x p_short_iip_backward,  long_iip = N x (1 - p_short_iip_backward)
#   averaged by calendar month (Jan-Dec). short/long IIP is a related-but-distinct
#   concept from same/cross-year in Fig 6, so the same colour tones are kept
#   (steel / terracotta).  "short IIP" here means I_H (infectee incubation period)
#   below the 70-day threshold used in compute_iip_probs.
#
#   Input : pipeline/processed/iip_probs_a0.3_b0.3.csv
#
# Run: Rscript figures/figS4_short_long_iip.R

suppressMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })

.args <- commandArgs(FALSE); .fa <- sub("^--file=", "", .args[grep("^--file=", .args)])
SCRIPT_DIR <- if (length(.fa)) dirname(normalizePath(.fa)) else getwd()
source(file.path(SCRIPT_DIR, "_common.R"))
OUT <- file.path(FIGDIR, "S4_Fig")

iip <- read_csv(file.path(DATA, "iip_probs_a0.3_b0.3.csv"), show_col_types = FALSE)
g <- iip %>% filter(!is.na(p_short_iip_backward)) %>%
  mutate(short = N * p_short_iip_backward, long = N * (1 - p_short_iip_backward)) %>%
  group_by(month) %>%
  summarise(short = mean(short), long = mean(long), .groups = "drop") %>%
  pivot_longer(c(short, long), names_to = "iip", values_to = "val") %>%
  mutate(iip = factor(iip, levels = c("short", "long"),
                      labels = c("short IIP", "long IIP")))
LV  <- c("short IIP", "long IIP")
pal <- setNames(c(STEEL, TERRACOTTA), LV)

p <- ggplot(g, aes(month, val, colour = iip, group = iip)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.8) +
  scale_colour_manual(values = pal) +
  scale_x_continuous("Month", breaks = 1:12, labels = month.abb) +
  scale_y_continuous("Mean weekly cases / month",
                     expand = expansion(mult = c(0.08, 0.05))) +  # bottom pad so y=0 points are not clipped
  theme_minimal_hgrid(font_size = 11, font_family = "Arial") +
  theme(axis.line.x = element_line(linewidth = 0.4),
        axis.title.x = element_text(margin = margin(t = 8)),
        axis.title.y = element_text(margin = margin(r = 8)),
        legend.position = "top", legend.title = element_blank(),
        legend.justification = "left", plot.margin = margin(6, 10, 6, 6))

save_fig(p, OUT, width = 5.8, height = 3.5)
