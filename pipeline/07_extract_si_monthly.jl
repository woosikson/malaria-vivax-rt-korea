#!/usr/bin/env julia
#
# Step 7 — monthly forward serial-interval distribution (input for Fig 2).
#   Group the per-t.4 forward SI distributions by calendar month of t.1 (= the
#   entry date, infector onset, restricted to 2013–2025). Within each month, sum
#   sample counts per si_week and normalise to a probability (sums to 1 per month).
#
#   Input  : pipeline/processed/forward_infector_si_a{REF}.jld2  (REF = α=β=0.3, from step 02)
#   Output : pipeline/processed/si_monthly_a{REF}.csv  (month,si_week,N,prob)
#            → read by the R figure scripts from pipeline/processed/si_monthly_a0.3_b0.3.csv
#
# Run: julia pipeline/07_extract_si_monthly.jl

using JLD2, DataFrames, Dates, CSV
include(joinpath(@__DIR__, "_config.jl"))

const T1_YEAR_MIN = 2013; const T1_YEAR_MAX = 2025

fn = joinpath(PROC, "forward_infector_si_a$(REF_A)_b$(REF_B).jld2")
si_all = jldopen(fn, "r"; iotype=IOStream) do file; file["si_dist"]; end
println("Loaded $(length(si_all)) entries from $(basename(fn))")

monthly_counts = Dict{Int, Dict{Int, Int}}()
for entry in si_all
    d = Date(entry.date); yr = year(d)
    (T1_YEAR_MIN <= yr <= T1_YEAR_MAX) || continue
    m = month(d)
    haskey(monthly_counts, m) || (monthly_counts[m] = Dict{Int, Int}())
    for row in eachrow(entry.si_dist)
        row.N > 0 && (monthly_counts[m][row.si_week] = get(monthly_counts[m], row.si_week, 0) + row.N)
    end
end

rows = NamedTuple[]
for m in 1:12
    (haskey(monthly_counts, m) && !isempty(monthly_counts[m])) || continue
    counts = monthly_counts[m]
    si_sorted = sort(collect(keys(counts)))
    ns = [counts[s] for s in si_sorted]; total = sum(ns)
    for k in eachindex(si_sorted)
        push!(rows, (month = m, si_week = si_sorted[k], N = ns[k], prob = ns[k] / total))
    end
end

out = joinpath(PROC, "si_monthly_a$(REF_A)_b$(REF_B).csv")
CSV.write(out, DataFrame(rows))
println("Saved: $(basename(out))  ($(length(rows)) rows)")
