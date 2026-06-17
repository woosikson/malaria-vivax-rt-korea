#!/usr/bin/env Rscript
#
# Fig 6 — Weekly cases I_t and R_t decomposed into same-year vs cross-year SI.
#   (A) case decomposition (backward SI, 2013-2025):
#         same  = round(N x ratio_backward),  cross = N - same
#   (B) R_t decomposition (forward SI, 2015-2023):
#         R_same = R_mean x ratio_forward,    R_cross = R_mean x (1 - ratio_forward)
#   same-year SI  = year(t0) == year(t4)  (infector & infectee onset in same year)
#   cross-year SI = different years (over-wintering / long latency crossing a year)
#
#   (A) cumulative area (histogram-like); (B) same/cross stacked area, broken in the
#   off-season (R_mean = 0) so winter gaps stay empty.
#   Colours: same-year = steel #3D6A93, cross-year = terracotta #C56B4A.  R=1 grey dashed.
#
#   Inputs : data/raw_inputs/malaria_kdca.csv
#            pipeline/processed/CI_R_case_a0.3_b0.3.csv
#            pipeline/processed/ratio_same_year_{backward,forward}_a0.3_b0.3.csv
#
# Run: Rscript figures/fig6_cases_and_Rt_decomposition.R

suppressMessages({ library(readr); library(dplyr); library(tidyr); library(lubridate); library(ggplot2) })

.args <- commandArgs(FALSE); .fa <- sub("^--file=", "", .args[grep("^--file=", .args)])
SCRIPT_DIR <- if (length(.fa)) dirname(normalizePath(.fa)) else getwd()
source(file.path(SCRIPT_DIR, "_common.R"))
OUT <- file.path(FIGDIR, "Fig6"); REF_COL <- "grey45"
LV  <- c("Same-year", "Cross-year")
pal <- setNames(c(STEEL, TERRACOTTA), LV)

ka <- read_csv(file.path(RAW,  "malaria_kdca.csv"),        show_col_types = FALSE)
rc <- read_csv(file.path(DATA, "CI_R_case_a0.3_b0.3.csv"), show_col_types = FALSE)
rb <- read_csv(file.path(DATA, "ratio_same_year_backward_a0.3_b0.3.csv"), show_col_types = FALSE)
rf <- read_csv(file.path(DATA, "ratio_same_year_forward_a0.3_b0.3.csv"),  show_col_types = FALSE)
ka$date <- as.Date(ka$date); rc$date <- as.Date(rc$date)
rb$date <- as.Date(rb$date); rf$date <- as.Date(rf$date)

XLIM <- as.Date(c("2013-01-01", "2026-01-01"))

# (A) case decomposition
ca <- ka %>% left_join(rb, by = "date") %>%
  filter(date >= XLIM[1], date < XLIM[2]) %>%
  mutate(ratio = coalesce(ratio, 0),
         `Same-year` = round(N * ratio), `Cross-year` = N - round(N * ratio)) %>%
  pivot_longer(c(`Same-year`, `Cross-year`), names_to = "comp", values_to = "val") %>%
  mutate(comp = factor(comp, levels = LV))

# (B) R_t decomposition over the full estimable window (2015-2023). Off-season
# R_mean = 0 -> both components 0 (flat, no fill); winter weeks with cases
# (R_mean > 0) appear as area spikes, matching panel (A) of Fig 3.
cb <- rc %>% select(date, R_mean) %>% left_join(rf, by = "date") %>%
  filter(date >= as.Date("2015-01-01"), date < as.Date("2024-01-01")) %>%
  arrange(date) %>%
  mutate(ratio = coalesce(ratio, 0),
         `Same-year` = R_mean * ratio, `Cross-year` = R_mean * (1 - ratio)) %>%
  pivot_longer(c(`Same-year`, `Cross-year`), names_to = "comp", values_to = "val") %>%
  mutate(comp = factor(comp, levels = LV))

xs <- function(title) scale_x_date(title, date_breaks = "1 year", date_labels = "%Y",
                                   limits = XLIM, expand = c(0.01, 0))
th_top <- theme_minimal_hgrid(font_size = 10, font_family = "Arial") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), plot.margin = margin(2, 8, 2, 4),
        axis.title.y = element_text(margin = margin(r = 8)),
        legend.position = "top", legend.title = element_blank(),
        legend.justification = "left")
th_bot <- theme_minimal_hgrid(font_size = 10, font_family = "Arial") +
  theme(axis.line.x = element_line(linewidth = 0.4), plot.margin = margin(2, 8, 4, 4),
        axis.title.x = element_text(margin = margin(t = 8)),
        axis.title.y = element_text(margin = margin(r = 8)),
        legend.position = "none")

RLAB <- expression(italic(R)[t])
p_a <- ggplot(ca, aes(date, val, fill = comp)) +
  geom_area(position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = pal) +
  xs(NULL) + scale_y_continuous("Weekly cases", expand = expansion(mult = c(0, 0.05))) +
  th_top
p_b <- ggplot(cb, aes(date, val, fill = comp)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = REF_COL, linewidth = 0.5) +
  geom_area(position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = pal) + xs("Date") +
  scale_y_continuous(RLAB, breaks = seq(0, 2, 0.5), expand = expansion(mult = c(0, 0.04))) +
  coord_cartesian(ylim = c(0, 2.05)) + th_bot   # stacked total (R_mean) max ~ 1.99 included

fig <- plot_grid(p_a, p_b, ncol = 1, align = "v", axis = "lr", rel_heights = c(1, 1),
                 labels = c("A", "B"), label_fontfamily = "Arial",
                 label_size = 13, label_fontface = "bold")
save_fig(fig, OUT, width = 7.5, height = 4.8)
