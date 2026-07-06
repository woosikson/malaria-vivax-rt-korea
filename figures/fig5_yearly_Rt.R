#!/usr/bin/env Rscript
#
# Fig 5 — Annual mean case reproduction number R_t (alpha=beta=0.3, 2015-2023).
#   Same scheme as Fig 4 / S2 Fig: black line + points, horizontal grid only,
#   grey dashed R=1 reference, size 5.8 x 3.5 in.
#   Aggregation: weekly R_t averaged per year over weeks with cases (N>0).
#
#   Inputs : data/raw_inputs/malaria_kdca.csv
#            pipeline/processed/CI_R_case_a0.3_b0.3.csv
#
# Run: Rscript figures/fig5_yearly_Rt.R

suppressMessages({ library(readr); library(dplyr); library(lubridate); library(ggplot2) })

.args <- commandArgs(FALSE); .fa <- sub("^--file=", "", .args[grep("^--file=", .args)])
SCRIPT_DIR <- if (length(.fa)) dirname(normalizePath(.fa)) else getwd()
source(file.path(SCRIPT_DIR, "_common.R"))
OUT <- file.path(FIGDIR, "Fig5"); REF_COL <- "grey50"

rc <- read_csv(file.path(DATA, "CI_R_case_a0.3_b0.3.csv"), show_col_types = FALSE)
ka <- read_csv(file.path(RAW,  "malaria_kdca.csv"),        show_col_types = FALSE)
rc$date <- as.Date(rc$date); ka$date <- as.Date(ka$date)
ann <- rc %>% left_join(ka, by = "date") %>%
  filter(date >= as.Date("2015-01-01"), date < as.Date("2024-01-01"), N > 0, R_mean > 0) %>%
  mutate(year = year(date)) %>% group_by(year) %>%
  summarise(R = mean(R_mean), .groups = "drop")

YLAB <- expression("Annual mean " * italic(R)[t])
p <- ggplot(ann, aes(year, R)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = REF_COL, linewidth = 0.5) +
  geom_line(colour = "black", linewidth = 0.7) +
  geom_point(colour = "black", size = 2.2) +
  scale_x_continuous("Year", breaks = 2015:2023) +
  scale_y_continuous(YLAB, breaks = seq(1.0, 1.6, 0.2)) +
  coord_cartesian(ylim = c(0.9, 1.6)) +
  theme_minimal_hgrid(font_size = 11, font_family = "Arial") +
  theme(axis.line.x = element_line(linewidth = 0.4), plot.margin = margin(6, 10, 6, 6),
        axis.title.x = element_text(margin = margin(t = 8)),
        axis.title.y = element_text(margin = margin(r = 8)))
save_fig(p, OUT, width = 5.8, height = 3.5)
