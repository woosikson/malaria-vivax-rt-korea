#!/usr/bin/env Rscript
#
# Fig 3 — Weekly cases I_t and the weekly case reproduction number R_t (alpha=beta=0.3).
#   (A) weekly cases (2013-2025)   (B) R_t line + 95% interval (2015-2023, ref R=1).
#   The full date range is kept to make clear that the first/last two years have
#   cases but no estimable R_t.  R_t is drawn as a continuous line+ribbon over the
#   whole estimable window (2015-2023) so sparse winter weeks are not broken out.
#
#   Visual encoding (neutral; colour is not used for meaning here):
#     cases   = grey needles (no fill)
#     R_t mean= black solid line
#     R_t = 1 = grey dashed reference
#     interval= light-grey ribbon
#
#   Inputs : data/raw_inputs/malaria_kdca.csv
#            pipeline/processed/CI_R_case_a0.3_b0.3.csv
#
# Run: Rscript figures/fig3_cases_and_Rt.R

suppressMessages({ library(readr); library(dplyr); library(lubridate); library(ggplot2) })

.args <- commandArgs(FALSE); .fa <- sub("^--file=", "", .args[grep("^--file=", .args)])
SCRIPT_DIR <- if (length(.fa)) dirname(normalizePath(.fa)) else getwd()
source(file.path(SCRIPT_DIR, "_common.R"))
OUT <- file.path(FIGDIR, "Fig3")

CASE_COL <- "grey35"; MEAN_COL <- "black"; RIBBON <- "grey85"; REF_COL <- "grey50"

ka <- read_csv(file.path(RAW,  "malaria_kdca.csv"),        show_col_types = FALSE)
rc <- read_csv(file.path(DATA, "CI_R_case_a0.3_b0.3.csv"), show_col_types = FALSE)
ka$date <- as.Date(ka$date); rc$date <- as.Date(rc$date)
cases <- ka %>% filter(date >= as.Date("2013-01-01"), date < as.Date("2026-01-01"))
# Continuous line+CI over the full estimable window (2015-2023): rc has weekly
# values in the off-season too, so not breaking keeps sparse winter weeks drawn
# as the same solid line + ribbon as the mid-season (matches the original).
rj <- rc %>% left_join(ka, by = "date") %>%
  filter(date >= as.Date("2015-01-01"), date < as.Date("2024-01-01")) %>%
  arrange(date)

XLIM <- as.Date(c("2013-01-01", "2026-01-01"))
RLAB <- expression(italic(R)[t])
xs <- function(title) scale_x_date(title, date_breaks = "1 year", date_labels = "%Y",
                                   limits = XLIM, expand = c(0.01, 0))

th_top <- theme_minimal_hgrid(font_size = 10, font_family = "Arial") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), plot.margin = margin(4, 8, 2, 4),
        axis.title.y = element_text(margin = margin(r = 8)))
th_bot <- theme_minimal_hgrid(font_size = 10, font_family = "Arial") +
  theme(axis.line.x = element_line(linewidth = 0.4), plot.margin = margin(2, 8, 4, 4),
        axis.title.x = element_text(margin = margin(t = 8)),
        axis.title.y = element_text(margin = margin(r = 8)))

p_a <- ggplot(cases, aes(date)) +
  geom_segment(aes(xend = date, y = 0, yend = N), colour = CASE_COL, linewidth = 0.3) +
  xs(NULL) + scale_y_continuous("Weekly cases", expand = expansion(mult = c(0, 0.05))) +
  th_top

p_b <- ggplot(rj, aes(date)) +
  geom_ribbon(aes(ymin = R_ll, ymax = R_ul), fill = RIBBON) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = REF_COL, linewidth = 0.5) +
  geom_line(aes(y = R_mean), colour = MEAN_COL, linewidth = 0.35) +
  xs("Date") + scale_y_continuous(RLAB, breaks = seq(0, 5, 1),
                                  expand = expansion(mult = c(0, 0.03))) +
  coord_cartesian(ylim = c(0, 5)) + th_bot   # interval max R_ul = 5.0 (2021-11-10) included

fig <- plot_grid(p_a, p_b, ncol = 1, align = "v", axis = "lr",
                 rel_heights = c(1, 1), labels = c("A", "B"),
                 label_fontfamily = "Arial", label_size = 13, label_fontface = "bold")
save_fig(fig, OUT, width = 7.5, height = 4.6)
