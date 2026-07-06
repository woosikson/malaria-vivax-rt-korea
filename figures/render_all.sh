#!/usr/bin/env bash
#
# render_all.sh — render all 10 manuscript figures from the figure inputs in
#   pipeline/processed/ (run pipeline/run_pipeline.sh first).
#   Each script writes <name>.eps and <name>.tiff to output/fig/.
#
# Usage: bash figures/render_all.sh
set -euo pipefail
cd "$(dirname "$0")"

for s in fig1_serial_interval_diagram.R \
         fig2_serial_interval_pdf.R \
         fig3_cases_and_Rt.R \
         fig4_monthly_Rt_expected_cases.R \
         fig5_yearly_Rt.R \
         fig6_cases_and_Rt_decomposition.R \
         figS1_IM_plus_TMH_pdf.R \
         figS2_yearly_Rt_sensitivity.R \
         figS3_monthly_Rt_sensitivity.R \
         figS4_short_long_iip.R; do
  echo "== $s =="
  Rscript "$s"
done
echo "== Done. Figures are in output/fig/ =="
