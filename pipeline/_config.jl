#
# _config.jl — shared paths, parameters, and helpers for the Julia pipeline.
#
# The pipeline regenerates every figure input from the three small CSVs in
# data/raw_inputs/.  Heavy intermediates (per-week Monte-Carlo samples, ~GBs)
# are written under pipeline/processed/ (git-ignored).  The small figure-input
# CSVs it produces are also written under pipeline/processed/, where the R
# figure scripts read them.
#
# `include`d by every pipeline/NN_*.jl script.

using Dates

const PIPE_DIR = @__DIR__
const ROOT     = normpath(joinpath(PIPE_DIR, ".."))
const RAW      = joinpath(ROOT, "data", "raw_inputs")
const PROC     = joinpath(PIPE_DIR, "processed")

# ── input files (small, non-confidential, shipped in data/raw_inputs/) ────────
const CASE_CSV     = joinpath(RAW, "malaria_kdca.csv")           # date,N  weekly domestic cases
const MOSQUITO_CSV = joinpath(RAW, "mosquito_index_2013_2025.csv") # year,week,mosquito_index
const TEMP_CSV     = joinpath(RAW, "temperature_daily.csv")      # date,temp daily 4-region mean

# ── model parameters (h_HM = 1 variant; α scaled by ALPHA_SCALE) ──────────────
const ALPHA_SCALE = 12.0          # anchor: median of seasonal mosquito peaks
const ALPHA_LIST  = [0.2, 0.3, 0.4]
const BETA_LIST   = [0.2, 0.3, 0.4]
const REF_A       = 0.3           # representative (α,β) for single-combo outputs
const REF_B       = 0.3

mkpath(PROC)

# ── KDCA epi-week helpers (week containing Jan 1 = week 1, Sun–Sat) ───────────
function kdca_week_start(yr::Int, wk::Int)::Date
    jan1 = Date(yr, 1, 1)
    wday_jan1 = dayofweek(jan1) % 7            # Julia Mon=1..Sun=7 → Sun=0
    week1_sunday = jan1 - Day(wday_jan1)
    return week1_sunday + Day((wk - 1) * 7)
end

function kdca_week_number(d::Date)::Int
    yr = year(d)
    jan1 = Date(yr, 1, 1)
    wday_jan1 = dayofweek(jan1) % 7
    week1_sunday = jan1 - Day(wday_jan1)
    return div(Dates.value(d - week1_sunday), 7) + 1
end
