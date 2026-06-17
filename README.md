# Reproducibility code and data — *Plasmodium vivax* malaria reproduction numbers (Republic of Korea)

This repository reproduces every figure in the manuscript from the original
input data. It contains the input data, the code that builds all intermediate
data, and the code that draws the figures.

> **Manuscript:** Son W-S, Lim A-Y, Hwang D-U, Nah K. *A time-varying risk
> assessment framework for* P. vivax *malaria transmission in temperate
> settings: A case study of the Republic of Korea.*

The figures estimate the **case reproduction number** *Rₜ* of vivax malaria from
weekly Korean surveillance data, using a mosquito-mediated serial-interval (SI)
model whose components are *T_HM* (human→mosquito transmission delay), *I_M*
(mosquito latent period), *T_MH* (mosquito→human transmission delay), and *I_H*
(infectee incubation period), so that the serial interval is
*t₅ − t₁ = T_HM + I_M + T_MH + I_H*.

## Reproducing the figures

Everything is regenerated from the three small input CSVs in `data/raw_inputs/`.
There are two steps: run the Julia pipeline to build the intermediate data, then
run the R scripts to draw the figures.

**Run all commands below from the repository root** — that is, the top-level
folder of this repository, the one that directly contains `data/`, `figures/`,
`pipeline/`, and this `README.md`. The commands use paths relative to that
folder (e.g. `pipeline/...`, `figures/...`), so they will not find their files if
run from inside a subdirectory. After cloning/downloading, `cd` into that folder
first, for example:

```bash
cd malaria-vivax-rt-korea     # the folder containing data/, figures/, pipeline/
```

```bash
# (still in the repository root)

# 1. build all intermediate data from data/raw_inputs/
julia --project=pipeline -e 'using Pkg; Pkg.instantiate()'   # first time only
bash pipeline/run_pipeline.sh

# 2. render the ten figures into output/fig/
bash figures/render_all.sh
```

> **The pipeline is heavy.** Generating the Monte-Carlo serial-interval samples
> for all nine (α, β) combinations takes **hours** of CPU time and writes on the
> order of **10+ GB** to `pipeline/processed/`. Use as many threads as you have
> cores (`--threads=auto`, already set in the driver). Each step skips work that
> is already complete, so the run can be interrupted and resumed.

The pipeline writes its outputs (including the small figure-input CSVs) to
`pipeline/processed/`; the R scripts read the figure inputs from there.

## Figure map (manuscript ↔ files)

| Manuscript | Output (`output/fig/`) | Script (`figures/`) | Content |
|------------|------------------------|---------------------|---------|
| **Fig 1**  | `Fig1.eps/.tiff`   | `fig1_serial_interval_diagram.R`       | Serial-interval decomposition schematic (conceptual) |
| **Fig 2**  | `Fig2.eps/.tiff`   | `fig2_serial_interval_pdf.R`           | Monthly forward serial-interval distribution (ridgeline) |
| **Fig 3**  | `Fig3.eps/.tiff`   | `fig3_cases_and_Rt.R`                  | Weekly cases (A) and weekly *Rₜ* with 95% interval (B) |
| **Fig 4**  | `Fig4.eps/.tiff`   | `fig4_yearly_Rt.R`                     | Annual mean *Rₜ* |
| **Fig 5**  | `Fig5.eps/.tiff`   | `fig5_monthly_Rt_expected_cases.R`     | Monthly mean *Rₜ* and expected secondary cases (dual axis) |
| **Fig 6**  | `Fig6.eps/.tiff`   | `fig6_cases_and_Rt_decomposition.R`    | Cases (A) and *Rₜ* (B) split into same-year / cross-year SI |
| **S1 Fig** | `S1_Fig.eps/.tiff` | `figS1_IM_plus_TMH_pdf.R`              | *I_M + T_MH* distribution, (α, β) sensitivity grid |
| **S2 Fig** | `S2_Fig.eps/.tiff` | `figS2_yearly_Rt_sensitivity.R`        | Annual mean *Rₜ*, (α, β) sensitivity grid |
| **S3 Fig** | `S3_Fig.eps/.tiff` | `figS3_monthly_Rt_sensitivity.R`       | Monthly mean *Rₜ*, (α, β) sensitivity grid |
| **S4 Fig** | `S4_Fig.eps/.tiff` | `figS4_short_long_iip.R`               | Monthly mean short-IIP / long-IIP cases |

(α, β) is the (human→mosquito, mosquito→human) transmission-hazard pair; α is
internally scaled by 12 (the median of the seasonal mosquito-density peaks).
The representative value used throughout the main text is α = β = 0.3; the
"sensitivity" figures (S1–S3) sweep {0.2, 0.3, 0.4}².

## Repository layout

```
.
├── data/
│   └── raw_inputs/        # the 3 small CSVs the whole pipeline runs from
│       ├── malaria_kdca.csv             # date, N — weekly domestic vivax cases
│       ├── mosquito_index_2013_2025.csv # year, week, mosquito_index — weekly mosquito surveillance
│       └── temperature_daily.csv        # date, temp — daily 4-region mean temperature
├── pipeline/             # Julia: raw_inputs → all intermediates (steps 01–08)
│   ├── processed/        # all pipeline outputs land here (git-ignored, can be large)
│   └── run_pipeline.sh
└── figures/              # R plotting code (one script per figure) + render_all.sh

output/fig/ is created by figures/render_all.sh (EPS + TIFF); it is not committed.
```

The figure-input CSVs produced under `pipeline/processed/`:

| File | Produced by | Used by |
|------|-------------|---------|
| `CI_R_case_a{α}_b{β}.csv` (9) | step 03 | Fig 3, 4, 5, 6, S2, S3 |
| `iip_probs_a0.3_b0.3.csv`     | step 04 | Fig 6(A), S4 |
| `si_monthly_a0.3_b0.3.csv`    | step 07 | Fig 2 |
| `ratio_same_year_{backward,forward}_a0.3_b0.3.csv` | step 06 | Fig 6 |
| `im_tmh_hist.csv`, `im_tmh_means.csv` | step 08 | S1 |

`malaria_kdca.csv` is also read directly by the plotting code (weekly case counts
for the Fig 3 / Fig 6 case panels).

## Data flow

```
raw_inputs (mosquito, temperature)
      │  01_generate_SI_{forward,backward}      (Monte-Carlo SI samples, ~GBs)
      ▼
SI_{forward,backward}_a{α}_b{β}/  (per-week sample files)
      │
      ├─ 02_aggregate_SI_dist_forward ─▶ forward_infector_si_a{α}_b{β}.jld2
      │        │  03_compute_R_case (+ malaria_kdca) ─▶ CI_R_case_a{α}_b{β}.csv      → Fig 3,4,5,6,S2,S3
      │        └─ 07_extract_si_monthly           ─▶ si_monthly_a0.3_b0.3.csv         → Fig 2
      ├─ 04_compute_iip_probs (+ malaria_kdca)    ─▶ iip_probs_a0.3_b0.3.csv          → Fig 6(A), S4
      ├─ 05_compute_hist_IM_TMH ─▶ hist_IM_TMH.jld2 ─ 08_extract_im_tmh_pdf ─▶ im_tmh_*.csv → S1
      └─ 06_extract_same_year_ratio              ─▶ ratio_same_year_*_a0.3_b0.3.csv   → Fig 6
```

## Reproducibility notes

- **Plotting is deterministic.** Given the same figure inputs, the R scripts
  reproduce the figures exactly.
- **The intermediate-data extractors are deterministic** (steps 04, 06, 07, 08):
  they produce identical CSVs from the serial-interval samples.
- **Serial-interval sampling (step 01) is seeded** per target week
  (`MersenneTwister(index)`), so it is reproducible given the same Julia and
  package versions.
- **Rₜ estimation (step 03) is stochastic** (10 000 Binomial Monte-Carlo trials).
  A fixed seed makes re-runs reproducible among themselves; over 10 000 trials
  the means are stable to ~3 decimals, leaving the figures visually identical.

## Software environment

- **R** ≥ 4.5 with `ggplot2` (≥ 4.0), `cowplot`, `ggridges`, `ragg`, `readr`,
  `dplyr`, `tidyr`, `lubridate`.
  ```r
  install.packages(c("ggplot2","cowplot","ggridges","ragg","readr","dplyr","tidyr","lubridate"))
  ```
  Fonts use **Arial**; on macOS install **XQuartz** so the cairo devices embed
  Arial as vector text in the EPS files.
- **Julia** ≥ 1.9; dependencies are pinned in `pipeline/Project.toml`
  (`Pkg.instantiate()` once). Run with `--threads=auto`.

## Data provenance and privacy

The three files in `data/raw_inputs/` are aggregate, non-identifying series:
weekly domestic case **counts** (KDCA), the weekly **mosquito surveillance**
index, and **daily mean temperature**. No patient-level line lists are included
or required — the entire pipeline runs from these aggregates.

## License

- **Code** (everything in `figures/` and `pipeline/`): MIT License — see `LICENSE`.
- **Data** (`data/raw_inputs/`): Creative Commons Attribution 4.0 (CC-BY 4.0).

Both run on freely available software only (R and Julia, with the packages listed
above); no proprietary or otherwise unobtainable software is required.
