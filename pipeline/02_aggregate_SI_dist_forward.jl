#!/usr/bin/env julia
#
# Step 2 (forward) — aggregate per-t.4 sample files into per-t.4 SI distributions.
#   For each (α,β) and each t.4 (Wednesday):
#     1. load pipeline/processed/SI_forward_a{α}_b{β}/{t4}.jld2
#     2. drop rows with missing t.0
#     3. keep rows with code t.1 and t.3 both in KDCA weeks 14–44
#     4. si_week = round((t.0 − t.4)/7); count → prob (normalised within the t.4)
#
#   Output : pipeline/processed/forward_infector_si_a{α}_b{β}.jld2
#            (a vector of (date, si_dist::DataFrame[si_week,N,prob]) tuples)
#   Used by : 03_compute_R_case.jl and 07_extract_si_monthly.jl
#
# Run: julia --threads=auto pipeline/02_aggregate_SI_dist_forward.jl

using DataFrames, Dates, JLD2
include(joinpath(@__DIR__, "_config.jl"))

println("Julia threads: $(Threads.nthreads())")

for alpha_orig in ALPHA_LIST, beta_val in BETA_LIST
    dn    = joinpath(PROC, "SI_forward_a$(alpha_orig)_b$(beta_val)")
    outfn = joinpath(PROC, "forward_infector_si_a$(alpha_orig)_b$(beta_val).jld2")
    if isfile(outfn); println("[α=$alpha_orig, β=$beta_val] SKIP (output exists)"); continue; end
    if !isdir(dn);    println("[α=$alpha_orig, β=$beta_val] SKIP (input dir missing)"); continue; end

    files      = sort(filter(f -> endswith(f, ".jld2"), readdir(dn)))
    file_paths = [joinpath(dn, f) for f in files]
    n_files    = length(file_paths)
    println("[α=$alpha_orig, β=$beta_val] Processing $n_files files...")
    t_start = time()
    results = Vector{Any}(nothing, n_files)

    empty_dist(ds) = (date = ds, si_dist = DataFrame(si_week = Int[0], N = Int[0], prob = Float64[0.0]))

    Threads.@threads for fi in 1:n_files
        fn = file_paths[fi]; date_str = splitext(basename(fn))[1]
        dt = jldopen(fn, "r"; iotype=IOStream) do file; file["dt"]; end
        if nrow(dt) == 0; results[fi] = empty_dist(date_str); continue; end
        valid = .!ismissing.(dt.t0)
        if !any(valid); results[fi] = empty_dist(date_str); continue; end
        dt_valid = dt[valid, :]
        t1_week = [kdca_week_number(d) for d in dt_valid.t1]
        t3_week = [kdca_week_number(d) for d in dt_valid.t3]
        week_mask = (14 .<= t3_week .<= 44) .& (14 .<= t1_week .<= 44)
        if !any(week_mask); results[fi] = empty_dist(date_str); continue; end
        dt_f = dt_valid[week_mask, :]
        si_weeks = [round(Int, Dates.value(dt_f.t0[i] - dt_f.t4[i]) / 7) for i in 1:nrow(dt_f)]
        count_dict = Dict{Int, Int}()
        for sw in si_weeks; count_dict[sw] = get(count_dict, sw, 0) + 1; end
        si_sorted = sort(collect(keys(count_dict)))
        counts    = [count_dict[s] for s in si_sorted]
        results[fi] = (date = date_str,
            si_dist = DataFrame(si_week = si_sorted, N = counts, prob = counts ./ sum(counts)))
    end

    jldopen(outfn, "w"; iotype=IOStream, compress=true) do file; file["si_dist"] = results; end
    println("[α=$alpha_orig, β=$beta_val] Done ($(round((time()-t_start)/60, digits=1)) min) -> $(basename(outfn))")
end
println("\nAll 9 forward SI distributions computed.")
