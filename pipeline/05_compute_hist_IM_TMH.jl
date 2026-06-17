#!/usr/bin/env julia
#
# Step 5 — I_M + T_MH distribution per (α,β) combo (9 combinations).
#   For each forward sample with t.1, t.2 non-missing and t.1 (= filename date)
#   in 2013–2025: T_MH = t.1 − t.2, value = I_M + T_MH, binned over 0..60 days.
#   Counts are normalised to a PDF (sums to 1); the mean is the exact sample mean.
#
#   Input  : pipeline/processed/SI_forward_a{α}_b{β}/   (per-week sample files)
#   Output : pipeline/processed/hist_IM_TMH.jld2  (counts, means, alpha_list, beta_list)
#   Used by: 08_extract_im_tmh_pdf.jl  (→ S1 Fig inputs)
#
# Run: julia --threads=auto pipeline/05_compute_hist_IM_TMH.jl

using DataFrames, Dates, JLD2
include(joinpath(@__DIR__, "_config.jl"))

println("Julia threads: $(Threads.nthreads())")
const T1_YEAR_MIN = 2013; const T1_YEAR_MAX = 2025

all_counts = Dict{String, Vector{Float64}}()
all_means  = Dict{String, Float64}()

for alpha_orig in ALPHA_LIST, beta_val in BETA_LIST
    dn = joinpath(PROC, "SI_forward_a$(alpha_orig)_b$(beta_val)")
    all_files = sort(filter(f -> endswith(f, ".jld2"), readdir(dn)))
    files = filter(all_files) do f
        d = Date(splitext(f)[1]); T1_YEAR_MIN <= year(d) <= T1_YEAR_MAX
    end
    file_paths = [joinpath(dn, f) for f in files]; n_files = length(file_paths)
    println("  [α=$alpha_orig, β=$beta_val] $n_files files")

    file_counts = [zeros(Float64, 61) for _ in 1:n_files]
    file_sums   = zeros(Float64, n_files); file_ns = zeros(Int, n_files)
    Threads.@threads for fi in 1:n_files
        dt = jldopen(file_paths[fi], "r"; iotype=IOStream) do file; file["dt"]; end
        nrow(dt) == 0 && continue
        for i in 1:nrow(dt)
            t1_i = dt.t1[i]; t2_i = dt.t2[i]
            (ismissing(t1_i) || ismissing(t2_i)) && continue
            IM_TMH = dt.I_M[i] + Dates.value(t1_i - t2_i)
            if 0 <= IM_TMH <= 60
                file_counts[fi][IM_TMH + 1] += 1.0
                file_sums[fi] += IM_TMH; file_ns[fi] += 1
            end
        end
    end

    counts = zeros(Float64, 61)
    for fi in 1:n_files; counts .+= file_counts[fi]; end
    total_count = sum(counts); total_count > 0 && (counts ./= total_count)
    total_n = sum(file_ns)
    key = "a$(alpha_orig)_b$(beta_val)"
    all_counts[key] = counts
    all_means[key]  = total_n > 0 ? sum(file_sums) / total_n : 0.0
end

outpath = joinpath(PROC, "hist_IM_TMH.jld2")
jldopen(outpath, "w"; iotype=IOStream) do file
    file["counts"] = all_counts; file["means"] = all_means
    file["alpha_list"] = ALPHA_LIST; file["beta_list"] = BETA_LIST
end
println("Saved: $(basename(outpath))")
