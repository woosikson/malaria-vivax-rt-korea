#!/usr/bin/env bash
#
# run_pipeline.sh — regenerate every figure input from data/raw_inputs/.
#
#   WARNING: steps 1–2 are heavy. Generating the Monte-Carlo serial-interval
#   samples for all 9 (alpha, beta) combinations takes HOURS of CPU time and
#   writes on the order of 10+ GB to pipeline/processed/. Use as many threads as
#   you have cores. You can interrupt and re-run: every step skips work that is
#   already complete.
#
#   The small figure-input CSVs land in pipeline/processed/, where the R figure
#   scripts (figures/render_all.sh) read them.
#
# Usage: bash pipeline/run_pipeline.sh
set -euo pipefail
cd "$(dirname "$0")"

JL=(julia --project=. --threads=auto)

echo "== Step 1: Monte-Carlo SI sampling (forward + backward) — HOURS, ~GBs =="
"${JL[@]}" 01_generate_SI_forward.jl
"${JL[@]}" 01_generate_SI_backward.jl

echo "== Step 2: aggregate forward SI distributions =="
"${JL[@]}" 02_aggregate_SI_dist_forward.jl

echo "== Step 3: weekly R_t with 95% interval (9 combos) =="
"${JL[@]}" 03_compute_R_case.jl

echo "== Step 4: short/long IIP probabilities (alpha=beta=0.3) =="
"${JL[@]}" 04_compute_iip_probs.jl

echo "== Step 5: I_M + T_MH histograms (9 combos) =="
"${JL[@]}" 05_compute_hist_IM_TMH.jl

echo "== Steps 6-8: extract remaining figure inputs =="
"${JL[@]}" 06_extract_same_year_ratio.jl
"${JL[@]}" 07_extract_si_monthly.jl
"${JL[@]}" 08_extract_im_tmh_pdf.jl

echo "== Pipeline complete. Figure inputs are in pipeline/processed/ =="
