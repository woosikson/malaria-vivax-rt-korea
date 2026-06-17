#!/usr/bin/env julia
#
# Step 4 — short/long IIP probability table for the domestic weekly case series.
#   IIP = I_H = infectee incubation period = (code t0 − code t1) in days.
#     short IIP : I_H <  IIP_THRESHOLD (70 d)      long IIP : I_H >= threshold
#   Sample filter (matches the SI decomposition): all 5 timepoints non-missing and
#   code t1, t3 both in KDCA weeks 14–44.
#   Backward : one file per target date (code t0 = target).
#   Forward  : scan all files; key each filtered sample by its own code t0.
#
#   Inputs : data/raw_inputs/malaria_kdca.csv
#            pipeline/processed/SI_backward_a{REF}/, SI_forward_a{REF}/  (REF = α=β=0.3)
#   Output : pipeline/processed/iip_probs_a{REF}.csv
#            → read by the R figure scripts from pipeline/processed/iip_probs_a0.3_b0.3.csv
#
# Run: julia --threads=auto pipeline/04_compute_iip_probs.jl

using CSV, DataFrames, Dates, JLD2, Base.Threads
include(joinpath(@__DIR__, "_config.jl"))

const A = REF_A; const B = REF_B; const IIP_THRESHOLD = 70
println("Julia threads: $(nthreads())  α=$A β=$B  IIP threshold=$IIP_THRESHOLD d")

df = DataFrame(CSV.File(CASE_CSV))
df.date = Date.(string.(df.date)); df.N = Int.(df.N)
df.week  = [kdca_week_number(d) for d in df.date]
df.month = [month(d) for d in df.date]
sort!(df, :date)
println("Loaded $(nrow(df)) weekly cases")

const target_dates    = df.date
const target_date_set = Set(target_dates)
const n_dates         = length(target_dates)

function valid_mask(dt::DataFrame)
    valid = .!ismissing.(dt.t0) .& .!ismissing.(dt.t1) .& .!ismissing.(dt.t2) .&
            .!ismissing.(dt.t3) .& .!ismissing.(dt.t4)
    any(valid) || return falses(nrow(dt))
    mask = copy(valid)
    @inbounds for i in eachindex(mask)
        mask[i] || continue
        if !(14 <= kdca_week_number(dt.t1[i]) <= 44 && 14 <= kdca_week_number(dt.t3[i]) <= 44)
            mask[i] = false
        end
    end
    return mask
end

# ── Backward: one file per target date ───────────────────────────────────────
const back_dn = joinpath(PROC, "SI_backward_a$(A)_b$(B)")
p_back = Vector{Union{Float64, Missing}}(undef, n_dates); n_back = zeros(Int, n_dates)
println("Backward: scanning $n_dates files in $(basename(back_dn))"); t_start = time()
@threads for i in 1:n_dates
    fn = joinpath(back_dn, "$(target_dates[i]).jld2")
    if !isfile(fn); p_back[i] = missing; continue; end
    dt = jldopen(fn, "r"; iotype=IOStream) do f; f["dt"]; end
    if nrow(dt) == 0; p_back[i] = missing; continue; end
    mask = valid_mask(dt); n_tot = count(mask)
    if n_tot == 0; p_back[i] = missing; continue; end
    t0_sel = dt.t0[mask]; t1_sel = dt.t1[mask]
    n_short = 0
    @inbounds for k in 1:n_tot
        Dates.value(t0_sel[k] - t1_sel[k]) < IIP_THRESHOLD && (n_short += 1)
    end
    p_back[i] = n_short / n_tot; n_back[i] = n_tot
end
println("  backward done in $(round((time()-t_start)/60, digits=2)) min")
df.p_short_iip_backward = p_back; df.n_back = n_back

# ── Forward: scan all files, aggregate by code t0 ────────────────────────────
const fwd_dn = joinpath(PROC, "SI_forward_a$(A)_b$(B)")
fwd_files = sort(filter(f -> endswith(f, ".jld2"), readdir(fwd_dn)))
const fwd_paths = [joinpath(fwd_dn, f) for f in fwd_files]
const n_files = length(fwd_paths)
println("Forward: scanning $n_files files in $(basename(fwd_dn))")
file_short = [Dict{Date, Int}() for _ in 1:n_files]
file_total = [Dict{Date, Int}() for _ in 1:n_files]
t_start = time()
@threads for fi in 1:n_files
    dt = jldopen(fwd_paths[fi], "r"; iotype=IOStream) do f; f["dt"]; end
    nrow(dt) == 0 && continue
    mask = valid_mask(dt); any(mask) || continue
    t0_sel = dt.t0[mask]; t1_sel = dt.t1[mask]
    ls = file_short[fi]; lt = file_total[fi]
    @inbounds for k in eachindex(t0_sel)
        tgt = t0_sel[k]; (tgt in target_date_set) || continue
        lt[tgt] = get(lt, tgt, 0) + 1
        Dates.value(tgt - t1_sel[k]) < IIP_THRESHOLD && (ls[tgt] = get(ls, tgt, 0) + 1)
    end
end
println("  forward scan done in $(round((time()-t_start)/60, digits=2)) min")

n_fwd_short = Dict{Date, Int}(d => 0 for d in target_dates)
n_fwd_total = Dict{Date, Int}(d => 0 for d in target_dates)
for fi in 1:n_files
    for (k, v) in file_total[fi]; n_fwd_total[k] += v; end
    for (k, v) in file_short[fi]; n_fwd_short[k] += v; end
end
p_fwd = Vector{Union{Float64, Missing}}(undef, n_dates); n_fwd = zeros(Int, n_dates)
for i in 1:n_dates
    d = target_dates[i]; n_tot = n_fwd_total[d]; n_fwd[i] = n_tot
    p_fwd[i] = n_tot > 0 ? n_fwd_short[d] / n_tot : missing
end
df.p_short_iip_forward = p_fwd; df.n_fwd = n_fwd

out_path = joinpath(PROC, "iip_probs_a$(A)_b$(B).csv")
CSV.write(out_path, df)
println("Saved: $(basename(out_path))")
