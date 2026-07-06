#!/usr/bin/env Rscript
#
# Fig 4 — Monthly (Jan-Dec) mean R_t and monthly mean expected secondary cases
#   (R_t x cases), on a dual axis (two lines).
#     left  axis : monthly mean R_t (line)
#     right axis : monthly mean (R_t x cases) = expected secondary cases / month
#       computed on weeks with cases (R != 0), then averaged by calendar month
#       over 2015-2023.
#   R_t = black, expected secondary cases = vermillion (#D55E00).  R=1 grey dashed.
#
#   Inputs : data/raw_inputs/malaria_kdca.csv
#            pipeline/processed/CI_R_case_a0.3_b0.3.csv
#
# Run: Rscript figures/fig4_monthly_Rt_expected_cases.R

suppressMessages({ library(readr); library(dplyr); library(lubridate); library(ggplot2) })

.args <- commandArgs(FALSE); .fa <- sub("^--file=", "", .args[grep("^--file=", .args)])
SCRIPT_DIR <- if (length(.fa)) dirname(normalizePath(.fa)) else getwd()
source(file.path(SCRIPT_DIR, "_common.R"))
OUT <- file.path(FIGDIR, "Fig4")

rc <- read_csv(file.path(DATA, "CI_R_case_a0.3_b0.3.csv"), show_col_types = FALSE)
ka <- read_csv(file.path(RAW,  "malaria_kdca.csv"),        show_col_types = FALSE)
rc$date <- as.Date(rc$date); ka$date <- as.Date(ka$date)
m <- rc %>% left_join(ka, by = "date") %>%
  filter(date >= as.Date("2015-01-01"), date < as.Date("2024-01-01"), N > 0, R_mean > 0) %>%
  mutate(month = month(date), prod = R_mean * N) %>%
  group_by(month) %>% summarise(R = mean(R_mean), sec = mean(prod), .groups = "drop")

# dual-axis mapping (right sec -> left R coordinate)
RLO <- 0.85; RHI <- 1.46; SECMAX <- 25
prim <- function(s) RLO + (s / SECMAX) * (RHI - RLO)
RLAB <- expression("Monthly mean " * italic(R)[t])
SLAB <- "Expected secondary cases / month"
xsc  <- scale_x_continuous("Month", breaks = 1:12, labels = month.abb,
                           expand = expansion(mult = c(0.03, 0.03)))
ysc  <- scale_y_continuous(RLAB, breaks = seq(0.9, 1.4, 0.1),
            sec.axis = sec_axis(~ (. - RLO) / (RHI - RLO) * SECMAX, name = SLAB,
                                breaks = seq(0, 25, 5)))
R_COL <- "black"; SEC_COL <- unname(OKABE["vermillion"])

p <- ggplot(m, aes(month)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  geom_line(aes(y = prim(sec)), colour = SEC_COL, linewidth = 0.7) +
  geom_point(aes(y = prim(sec)), colour = SEC_COL, size = 2.2, shape = 15) +
  geom_line(aes(y = R), colour = R_COL, linewidth = 0.7) +
  geom_point(aes(y = R), colour = R_COL, size = 2.2) +
  xsc + ysc + coord_cartesian(ylim = c(RLO, RHI)) +
  theme_minimal_hgrid(font_size = 11, font_family = "Arial") +
  theme(axis.line.x = element_line(linewidth = 0.4),
        axis.title.x = element_text(margin = margin(t = 8)),
        axis.title.y.left  = element_text(colour = R_COL, margin = margin(r = 8)),
        axis.text.y.left   = element_text(colour = R_COL),
        axis.title.y.right = element_text(colour = SEC_COL, margin = margin(l = 8)),
        axis.text.y.right  = element_text(colour = SEC_COL),
        plot.margin = margin(6, 10, 6, 6))

save_fig(p, OUT, width = 5.8, height = 3.5)
