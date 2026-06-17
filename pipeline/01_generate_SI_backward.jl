#!/usr/bin/env julia
#
# Step 1 (backward) — Monte-Carlo serial-interval sampling, infectee-anchored.
#   Backward flow from an infectee onset t.0 (= code t0 = concept t5):
#     t0 --I_H--> t1 (M→H bite) --T_MH--> t2 (mosquito infectious)
#        --I_M(backward)--> t3 (H→M bite) --T_HM--> t4 (infector onset, concept t1)
#   h_HM = 1; α ∈ {0.2,0.3,0.4}/ALPHA_SCALE, β ∈ {0.2,0.3,0.4} → 9 combos.
#
#   Inputs : data/raw_inputs/mosquito_index_2013_2025.csv, temperature_daily.csv
#   Output : pipeline/processed/SI_backward_a{α}_b{β}/{t0}.jld2  (one file per Wednesday t0)
#
#   Deterministic per-t0 seed. WARNING: long runtime (hours), ~1 GB per combo.
#   Only α=β=0.3 is strictly required downstream (same-year ratio + IIP probabilities),
#   but all 9 combos are generated to mirror the published methodology.
#
# Run: julia --threads=auto pipeline/01_generate_SI_backward.jl

using CSV, DataFrames, Dates, Distributions, Random, JLD2, StatsBase, SpecialFunctions
include(joinpath(@__DIR__, "_config.jl"))

println("Julia threads: $(Threads.nthreads())")

# ─── Load data ────────────────────────────────────────────────────────────────
println("Loading data...")
mos_weekly = CSV.read(MOSQUITO_CSV, DataFrame)
temp_data  = CSV.read(TEMP_CSV, DataFrame)
temp_data.date = Date.(string.(temp_data.date))

mos_N_dict = Dict{Date, Float64}()
for r in eachrow(mos_weekly)
    yr  = Int(r.year); wk = Int(r.week); val = Float64(r.mosquito_index)
    start_d = kdca_week_start(yr, wk)
    for off in 0:6
        mos_N_dict[start_d + Day(off)] = val
    end
end
const mos_N    = mos_N_dict
const temp_val = Dict{Date, Union{Float64, Missing}}(r.date => r.temp for r in eachrow(temp_data))

get_mos_N(d::Date)::Float64 = get(mos_N, d, 0.0)
function get_dev_rate(d::Date)::Float64
    t = get(temp_val, d, missing)
    (ismissing(t) || isnan(t) || t < 14.5) && return 0.0
    return (t - 14.5) / 105.0
end

# ─── Forward & backward I_M lookups ──────────────────────────────────────────
const all_dates       = collect(Date(2013,1,1):Day(1):Date(2025,12,31))
const n_all_dates     = length(all_dates)
const temp_end_date   = maximum(temp_data.date)
const temp_start_date = minimum(temp_data.date)

println("Computing forward I_M lookup...")
im_fwd_results = Vector{Union{Date, Nothing}}(nothing, n_all_dates)
Threads.@threads for idx in 1:n_all_dates
    d = all_dates[idx]; cum = 0.0; dd = d
    while dd <= temp_end_date
        cum += get_dev_rate(dd); cum > 1.0 && (im_fwd_results[idx] = dd; break); dd += Day(1)
    end
end
const im_fwd = Dict{Date, Union{Date, Nothing}}(all_dates[i] => im_fwd_results[i] for i in 1:n_all_dates)
println("  Forward I_M done. NA: $(count(isnothing, im_fwd_results))")

println("Computing backward I_M lookup...")
im_back_results = Vector{Union{Date, Nothing}}(nothing, n_all_dates)
Threads.@threads for idx in 1:n_all_dates
    d = all_dates[idx]; cum = 0.0; dd = d
    while dd >= temp_start_date
        cum += get_dev_rate(dd); cum > 1.0 && (im_back_results[idx] = dd; break); dd -= Day(1)
    end
end
const im_back = Dict{Date, Union{Date, Nothing}}(all_dates[i] => im_back_results[i] for i in 1:n_all_dates)
println("  Backward I_M done. NA: $(count(isnothing, im_back_results))")

# ─── t.3 candidates and possible t.1 ─────────────────────────────────────────
list_t3 = Date[]
for yr in 2013:2025
    start_d = kdca_week_start(yr, 14); end_d = kdca_week_start(yr, 44) + Day(6)
    append!(list_t3, collect(start_d:Day(1):end_d))
end
println("Computing possible_t.1...")
possible_t1_set = Set{Date}()
for u in list_t3
    t2 = get(im_fwd, u, nothing)
    isnothing(t2) && continue
    year(u) != year(t2) && continue
    I_M = Dates.value(t2 - u); T_MH_max = 60 - I_M
    T_MH_max < 0 && continue
    for j in 0:T_MH_max
        t1_try = t2 + Day(j)
        14 <= kdca_week_number(t1_try) <= 44 && push!(possible_t1_set, t1_try)
    end
end
println("  possible_t.1: $(length(possible_t1_set)) dates")

# ─── Parameters ──────────────────────────────────────────────────────────────
const t0_list = collect(Date(2013,1,2):Day(7):Date(2025,12,31))
println("t.0 list: $(length(t0_list)) dates")
const n_I_H = 500; const n_T_MH = 30; const n_T_HM = 100
const p_mix = 0.7423; const k_gamma = 22.8197; const theta_gamma = 1.11405
const mu_ln = 5.78509; const sigma_ln = 0.140988

# ─── Core backward function ──────────────────────────────────────────────────
function generate_SI_backward_for_t0(t0::Date, alpha::Float64, beta::Float64, rng::AbstractRNG)
    I_H_vals = 1:500
    Pr_I_H = Float64[
        p_mix     * v^(k_gamma-1) * exp(-v/theta_gamma) / (theta_gamma^k_gamma * gamma(k_gamma)) +
        (1-p_mix) * exp(-(log(v) - mu_ln)^2 / (2*sigma_ln^2)) / (v * sigma_ln * sqrt(2*pi))
        for v in I_H_vals]
    for j in 1:500
        (t0 - Day(j)) in possible_t1_set || (Pr_I_H[j] = 0.0)
    end
    sum(Pr_I_H) == 0.0 && return DataFrame()

    w_IH = Weights(Pr_I_H)
    I_H_samples = [sample(rng, I_H_vals, w_IH) for _ in 1:n_I_H]
    t1_samples  = [t0 - Day(v) for v in I_H_samples]

    n_012  = n_I_H * n_T_MH
    t0_012 = fill(t0, n_012)
    t1_012 = repeat(t1_samples, inner=n_T_MH)
    t2_012 = Vector{Union{Date, Missing}}(missing, n_012)

    for i in 1:n_I_H
        u = t1_samples[i]; base_idx = (i - 1) * n_T_MH
        N_TMH   = Float64[get_mos_N(u + Day(j)) for j in 0:60]
        Pr_T_MH = [(N_TMH[j+1] != 0.0 ? 1.0 : 0.0) * beta * exp(-beta * j) for j in 0:60]
        for j in 1:61
            Pr_T_MH[j] == 0.0 && continue
            T_MH_val = j - 1; t2_try = u - Day(T_MH_val)
            if t2_try < Date(2013,1,1) || t2_try > Date(2025,12,31); Pr_T_MH[j] = 0.0; continue; end
            t3_try = get(im_back, t2_try, nothing)
            if isnothing(t3_try); Pr_T_MH[j] = 0.0; continue; end
            I_M_val = Dates.value(t2_try - t3_try)
            if I_M_val + T_MH_val > 60; Pr_T_MH[j] = 0.0; continue; end
            year(t3_try) != year(t2_try) && (Pr_T_MH[j] = 0.0)
        end
        sum(Pr_T_MH) == 0.0 && continue
        w_TMH = Weights(Pr_T_MH)
        for k in 1:n_T_MH
            t2_012[base_idx + k] = u - Day(sample(rng, 0:60, w_TMH))
        end
    end

    valid_mask = .!ismissing.(t2_012)
    valid_idx  = findall(valid_mask)
    isempty(valid_idx) && return DataFrame()

    t3_012 = Vector{Union{Date, Missing}}(missing, n_012)
    for i in valid_idx
        t3_val = get(im_back, t2_012[i]::Date, nothing)
        isnothing(t3_val) ? (valid_mask[i] = false) : (t3_012[i] = t3_val)
    end
    valid_idx = findall(valid_mask .& .!ismissing.(t3_012))
    isempty(valid_idx) && return DataFrame()

    t0_v = t0_012[valid_idx]; t1_v = t1_012[valid_idx]
    t2_v = Date[t2_012[i]::Date for i in valid_idx]
    t3_v = Date[t3_012[i]::Date for i in valid_idx]

    unique_t3 = sort(unique(t3_v))
    t3_counts = Dict{Date, Int}()
    for d in t3_v; t3_counts[d] = get(t3_counts, d, 0) + 1; end

    t4_by_t3 = Dict{Date, Vector{Union{Date, Missing}}}()
    for u in unique_t3
        total = t3_counts[u] * n_T_HM
        N_vec = Float64[get_mos_N(u + Day(j)) for j in 0:365]
        cum_N = cumsum(N_vec)
        Pr_T_HM = alpha .* N_vec .* exp.(-alpha .* cum_N)
        if sum(Pr_T_HM) != 0.0
            w_THM = Weights(Pr_T_HM)
            T_HM_vals = [sample(rng, 0:365, w_THM) for _ in 1:total]
            t4_by_t3[u] = [u - Day(v) for v in T_HM_vals]
        else
            t4_by_t3[u] = fill(missing, total)
        end
    end

    sort_perm = sortperm(t3_v)
    t0_sorted = t0_v[sort_perm]; t1_sorted = t1_v[sort_perm]
    t2_sorted = t2_v[sort_perm]; t3_sorted = t3_v[sort_perm]

    t0_final = repeat(t0_sorted, inner=n_T_HM)
    t1_final = repeat(t1_sorted, inner=n_T_HM)
    t2_final = repeat(t2_sorted, inner=n_T_HM)
    t3_final = repeat(t3_sorted, inner=n_T_HM)
    t4_final = Vector{Union{Date, Missing}}(missing, length(t0_final))
    t3_pos   = Dict{Date, Int}(u => 0 for u in unique_t3)
    for i in eachindex(t3_final)
        u = t3_final[i]; pos = t3_pos[u] + 1; t3_pos[u] = pos
        t4_final[i] = t4_by_t3[u][pos]
    end
    I_M_final = [Dates.value(t2_final[i] - t3_final[i]) for i in eachindex(t2_final)]

    return DataFrame(t0 = t0_final, t1 = t1_final, t2 = t2_final,
                     t3 = t3_final, t4 = t4_final, I_M = I_M_final)
end

# ─── Run 9 parameter combinations ────────────────────────────────────────────
for alpha_orig in ALPHA_LIST, beta_val in BETA_LIST
    alpha = alpha_orig / ALPHA_SCALE
    dn = joinpath(PROC, "SI_backward_a$(alpha_orig)_b$(beta_val)")
    mkpath(dn)
    if count(f -> endswith(f, ".jld2"), readdir(dn)) == length(t0_list)
        println("[α=$alpha_orig, β=$beta_val] SKIP (complete)"); continue
    end
    println("[α=$alpha_orig, β=$beta_val] Starting ($(length(t0_list)) dates)...")
    t_start = time(); save_lock = ReentrantLock()
    Threads.@threads for ind in 1:length(t0_list)
        t0 = t0_list[ind]; fn = joinpath(dn, "$(t0).jld2")
        isfile(fn) && continue
        dt = generate_SI_backward_for_t0(t0, alpha, beta_val, MersenneTwister(ind))
        lock(save_lock) do
            jldopen(fn, "w"; iotype=IOStream, compress=true) do file; file["dt"] = dt; end
        end
    end
    println("[α=$alpha_orig, β=$beta_val] Done ($(round((time()-t_start)/60, digits=1)) min)")
end
println("\nAll 9 backward combinations completed.")
