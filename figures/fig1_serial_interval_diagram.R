#!/usr/bin/env Rscript
#
# Fig 1 — Serial-interval decomposition schematic for vector-borne transmission.
#   serial interval = t5 - t1 = T_HM + I_M + T_MH + I_H.
#   Host phase is colour-coded (Okabe-Ito):
#     in human  (t1->t2, t4->t5) : blue       #0072B2
#     in mosquito (t2->t3->t4)   : vermillion #D55E00
#   so the pathogen's H -> M -> H route is readable at a glance.
#   Conceptual diagram (no data input).
#
# Run: Rscript figures/fig1_serial_interval_diagram.R

suppressMessages({ library(ggplot2) })

.args <- commandArgs(FALSE); .fa <- sub("^--file=", "", .args[grep("^--file=", .args)])
SCRIPT_DIR <- if (length(.fa)) dirname(normalizePath(.fa)) else getwd()
source(file.path(SCRIPT_DIR, "_common.R"))            # OKABE, save_fig(), FIGDIR
OUT <- file.path(FIGDIR, "Fig1")

HUMAN <- unname(OKABE["blue"])        # in human
MOSQ  <- unname(OKABE["vermillion"])  # in mosquito
INK   <- "grey15"                     # event arrows, labels, axis

# ── coordinates ──────────────────────────────────────────────────────────────
X       <- 1:5                  # t1..t5 (evenly spaced; schematic, not to scale)
H_ARC   <- 0.42                 # arc height
ARR_TOP <- 0.92                 # event-arrow top

ev <- data.frame(
  x   = X,
  tl  = c("italic(t)[1]", "italic(t)[2]", "italic(t)[3]", "italic(t)[4]", "italic(t)[5]"),
  lab = c("Infector's\nsymptom onset", "Transmission\n(H → M)",
          "Mosquito\nbeing infectious", "Transmission\n(M → H)",
          "Infectee's\nsymptom onset")
)

# interval arcs (dashed) -------------------------------------------------------
arc_xy <- function(a, b, h, col, lab, n = 120) {
  t <- seq(0, pi, length.out = n)
  data.frame(x = (a + b)/2 - (b - a)/2 * cos(t), y = h * sin(t), grp = lab, col = col)
}
arc_spec <- list(
  list(1, 2, HUMAN, "T[HM]"), list(2, 3, MOSQ, "I[M]"),
  list(3, 4, MOSQ, "T[MH]"), list(4, 5, HUMAN, "I[H]")
)
arcs <- do.call(rbind, lapply(arc_spec, function(s) arc_xy(s[[1]], s[[2]], H_ARC, s[[3]], s[[4]])))
arc_lab <- data.frame(
  x   = sapply(arc_spec, function(s) (s[[1]] + s[[2]])/2),
  y   = H_ARC + 0.07,
  lab = sapply(arc_spec, `[[`, 4),
  col = sapply(arc_spec, `[[`, 3)
)

# host-phase brackets (below axis) --------------------------------------------
yb <- -0.40
spans <- data.frame(
  x0  = c(1, 2, 4), x1 = c(2, 4, 5),
  lab = c("in human", "in mosquito", "in human"),
  col = c(HUMAN, MOSQ, HUMAN)
)

# ── plot ─────────────────────────────────────────────────────────────────────
p <- ggplot() +
  geom_segment(aes(x = 0.45, xend = 5.55, y = 0, yend = 0), linewidth = 0.5, colour = INK) +
  geom_path(data = arcs, aes(x, y, group = grp, colour = col),
            linetype = "22", linewidth = 0.6, lineend = "round") +
  geom_text(data = arc_lab, aes(x, y, label = lab, colour = col),
            parse = TRUE, family = "Arial", size = 3.5, vjust = 0) +
  geom_segment(data = ev, aes(x = x, xend = x, y = ARR_TOP, yend = 0.015),
               colour = INK, linewidth = 0.45,
               arrow = arrow(length = unit(5, "pt"), type = "closed")) +
  geom_text(data = ev, aes(x, ARR_TOP + 0.07, label = lab),
            family = "Arial", size = 3.05, colour = INK, lineheight = 0.95, vjust = 0) +
  geom_text(data = ev, aes(x, -0.14, label = tl), parse = TRUE,
            family = "Arial", size = 3.9, colour = INK, vjust = 1) +
  geom_segment(data = spans, aes(x = x0, xend = x1, y = yb, yend = yb, colour = col), linewidth = 0.5) +
  geom_segment(data = spans, aes(x = x0, xend = x0, y = yb, yend = yb + 0.05, colour = col), linewidth = 0.5) +
  geom_segment(data = spans, aes(x = x1, xend = x1, y = yb, yend = yb + 0.05, colour = col), linewidth = 0.5) +
  geom_text(data = spans, aes((x0 + x1)/2, yb - 0.09, label = lab, colour = col),
            family = "Arial", size = 2.9, fontface = "italic", vjust = 1) +
  scale_colour_identity() +
  coord_cartesian(xlim = c(0.4, 5.6), ylim = c(-0.62, 1.42), clip = "off") +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(4, 8, 4, 8))

save_fig(p, OUT, width = 7.5, height = 2.7)
