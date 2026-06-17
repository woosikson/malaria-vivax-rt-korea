#!/usr/bin/env Rscript
#
# S3 Fig — Sensitivity analysis of monthly mean R_t over the (alpha, beta) grid (3x3).
#   9 combinations (alpha, beta) in {0.2, 0.3, 0.4}^2; monthly mean R_t
#   (Jan-Dec, 2015-2023, weeks with cases). The monthly pattern is nearly
#   invariant to alpha, beta (all combos 0.91-1.41, same shape) -> robust.
#   Same scheme as S2 Fig: facet_grid (beta-rows x alpha-cols), black line,
#   horizontal grid only, sub-panel aspect ratio = Fig 5 (5.8:3.5).
#
#   Inputs : data/raw_inputs/malaria_kdca.csv
#            pipeline/processed/CI_R_case_a{alpha}_b{beta}.csv  (9 files)
#
# Run: Rscript figures/figS3_monthly_Rt_sensitivity.R

suppressMessages({ library(readr); library(dplyr); library(lubridate); library(ggplot2) })

.args <- commandArgs(FALSE); .fa <- sub("^--file=", "", .args[grep("^--file=", .args)])
SCRIPT_DIR <- if (length(.fa)) dirname(normalizePath(.fa)) else getwd()
source(file.path(SCRIPT_DIR, "_common.R"))
OUT <- file.path(FIGDIR, "S3_Fig"); REF_COL <- "grey50"

ka <- read_csv(file.path(RAW, "malaria_kdca.csv"), show_col_types = FALSE)
ka$date <- as.Date(ka$date)
d <- bind_rows(lapply(c(0.2, 0.3, 0.4), function(a) bind_rows(lapply(c(0.2, 0.3, 0.4), function(b) {
  rc <- read_csv(file.path(DATA, sprintf("CI_R_case_a%s_b%s.csv", a, b)), show_col_types = FALSE)
  rc$date <- as.Date(rc$date)
  rc %>% left_join(ka, by = "date") %>%
    filter(date >= as.Date("2015-01-01"), date < as.Date("2024-01-01"), N > 0, R_mean > 0) %>%
    mutate(month = month(date)) %>% group_by(month) %>%
    summarise(R = mean(R_mean), .groups = "drop") %>% mutate(alpha = a, beta = b)
}))))

MAB  <- month.abb
YLAB <- expression("Monthly mean " * italic(R)[t])
xsc <- scale_x_continuous("Month", breaks = c(1, 4, 7, 10), labels = MAB[c(1, 4, 7, 10)])
ysc <- scale_y_continuous(YLAB, breaks = seq(0.9, 1.4, 0.1))
ref <- geom_hline(yintercept = 1, linetype = "dashed", colour = REF_COL, linewidth = 0.4)

p <- ggplot(d, aes(month, R)) + ref +
  geom_line(colour = "black", linewidth = 0.6) + geom_point(colour = "black", size = 1) +
  facet_grid(beta ~ alpha,
             labeller = label_bquote(cols = alpha == .(alpha), rows = beta == .(beta))) +
  xsc + ysc + coord_cartesian(ylim = c(0.88, 1.45)) +
  theme_minimal_hgrid(font_size = 10, font_family = "Arial") +
  theme(panel.spacing = unit(8, "pt"), axis.text = element_text(size = 8),
        axis.line.x = element_line(linewidth = 0.4),
        axis.title.x = element_text(margin = margin(t = 8)),
        axis.title.y = element_text(margin = margin(r = 8)),
        aspect.ratio = 3.5 / 5.8)
save_fig(p, OUT, width = 7.4, height = 5.2)
