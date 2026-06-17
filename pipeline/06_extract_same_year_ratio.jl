#!/usr/bin/env julia
#
# Step 6 — same-year SI ratio per week (input for the Fig 6 decomposition).
#   For each week's raw samples, fraction with year(t0) == year(t4) (same-year
#   transmission) among samples passing the analysis filter (code t1, t3 in KDCA
#   weeks 14–44):
#     backward (panel A): require all of t0..t4 non-missing
#     forward  (panel B): require t0 non-missing
#
#   Input  : pipeline/processed/SI_backward_a{REF}/, SI_forward_a{REF}/ (REF = α=β=0.3)
#   Output : pipeline/processed/ratio_same_year_{backward,forward}_a{REF}.csv  (date,ratio)
#            → read by the R figure scripts from pipeline/processed/ratio_same_year_{backward,forward}_a0.3_b0.3.csv
#
# Run: julia -t auto pipeline/06_extract_same_year_ratio.jl

using JLD2, DataFrames, Dates, CSV, Base.Threads
include(joinpath(@__DIR__, "_config.jl"))

# require_full = true → backward filter (all t0..t4); false → forward filter (t0 only)
function ratio_dir(dn::String, require_full::Bool)
    files = sort(filter(f -> endswith(f, ".jld2"), readdir(dn)))
    paths = [joinpath(dn, f) for f in files]; n = length(paths)
    out = Vector{Tuple{Date, Float64}}(undef, n)
    @threads for i in 1:n
        ds = splitext(basename(paths[i]))[1]
        dt = jldopen(paths[i], "r"; iotype=IOStream) do f; f["dt"]; end
        if nrow(dt) == 0; out[i] = (Date(ds), 0.0); continue; end
        valid = require_full ?
            (.!ismissing.(dt.t0) .& .!ismissing.(dt.t1) .& .!ismissing.(dt.t2) .&
             .!ismissing.(dt.t3) .& .!ismissing.(dt.t4)) :
            .!ismissing.(dt.t0)
        if !any(valid); out[i] = (Date(ds), 0.0); continue; end
        dv = dt[valid, :]
        t1w = [kdca_week_number(d) for d in dv.t1]
        t3w = [kdca_week_number(d) for d in dv.t3]
        wm = (14 .<= t3w .<= 44) .& (14 .<= t1w .<= 44)
        if !any(wm); out[i] = (Date(ds), 0.0); continue; end
        df = dv[wm, :]; nt = nrow(df); ns = 0
        for j in 1:nt
            year(df.t0[j]) == year(df.t4[j]) && (ns += 1)
        end
        out[i] = (Date(ds), ns / nt)
    end
    sort!(out, by = x -> x[1])
    DataFrame(date = [x[1] for x in out], ratio = [x[2] for x in out])
end

println("threads=$(nthreads())  computing backward ratio...")
@time back = ratio_dir(joinpath(PROC, "SI_backward_a$(REF_A)_b$(REF_B)"), true)
CSV.write(joinpath(PROC, "ratio_same_year_backward_a$(REF_A)_b$(REF_B).csv"), back)
println("computing forward ratio...")
@time fwd = ratio_dir(joinpath(PROC, "SI_forward_a$(REF_A)_b$(REF_B)"), false)
CSV.write(joinpath(PROC, "ratio_same_year_forward_a$(REF_A)_b$(REF_B).csv"), fwd)
println("done: $(nrow(back)) backward rows, $(nrow(fwd)) forward rows")
