#!/usr/bin/env julia
#
# Step 1 (forward) — Monte-Carlo serial-interval sampling, infector-anchored.
#   Forward flow from an infector onset t.1 (= code t4):
#     t4 --T_HM--> t3 (H→M bite) --I_M--> t2 (mosquito infectious)
#        --T_MH--> t1 (M→H bite) --I_H--> t0 (infectee onset, code t0 = concept t5)
#   h_HM = 1 (h_HM(s,t) = α·M(t)); α ∈ {0.2,0.3,0.4}/ALPHA_SCALE, β ∈ {0.2,0.3,0.4} → 9 combos.
#
#   Inputs : data/raw_inputs/mosquito_index_2013_2025.csv, temperature_daily.csv
#   Output : pipeline/processed/SI_forward_a{α}_b{β}/{t4}.jld2  (one file per Wednesday t4)
#
#   Deterministic: each t4 uses MersenneTwister(index), so runs are reproducible.
#   WARNING: long runtime (hours) and large output (~1 GB per combo).
#
# Run: julia --threads=auto pipeline/01_generate_SI_forward.jl

using CSV, DataFrames, Dates, Distributions, Random, JLD2, StatsBase
include(joinpath(@__DIR__, "_config.jl"))

println("Julia threads: $(Threads.nthreads())")

# ─── Load data ────────────────────────────────────────────────────────────────
println("Loading data...")
mos_weekly = CSV.read(MOSQUITO_CSV, DataFrame)
temp_data  = CSV.read(TEMP_CSV, DataFrame)
temp_data.date = Date.(string.(temp_data.date))

# ─── Daily mosquito index (weekly value broadcast to its 7 days) ──────────────
mos_N_dict = Dict{Date, Float64}()
for r in eachrow(mos_weekly)
    yr  = Int(r.year); wk = Int(r.week); val = Float64(r.mosquito_index)
    start_d = kdca_week_start(yr, wk)
    for off in 0:6
        mos_N_dict[start_d + Day(off)] = val
    end
end
const mos_N = mos_N_dict
println("mos_N: $(length(mos_N)) daily entries")

const temp_val = Dict{Date, Union{Float64, Missing}}(r.date => r.temp for r in eachrow(temp_data))

get_mos_N(d::Date)::Float64 = get(mos_N, d, 0.0)
function get_dev_rate(d::Date)::Float64
    t = get(temp_val, d, missing)
    (ismissing(t) || isnan(t) || t < 14.5) && return 0.0
    return (t - 14.5) / 105.0
end

# ─── Precompute forward I_M lookup (mosquito latent period via degree-days) ───
println("Computing I_M lookup table...")
const all_dates     = collect(Date(2013,1,1):Day(1):Date(2025,12,31))
const n_all_dates   = length(all_dates)
const temp_end_date = maximum(temp_data.date)
im_results = Vector{Union{Date, Nothing}}(nothing, n_all_dates)
Threads.@threads for idx in 1:n_all_dates
    d = all_dates[idx]; cum = 0.0; dd = d
    while dd <= temp_end_date
        cum += get_dev_rate(dd)
        if cum > 1.0
            im_results[idx] = dd; break
        end
        dd += Day(1)
    end
end
const im_lookup = Dict{Date, Union{Date, Nothing}}(all_dates[i] => im_results[i] for i in 1:n_all_dates)
println("I_M lookup done. NA count: $(count(isnothing, im_results))")

# ─── t.3 candidates (KDCA weeks 14–44) and possible t.4 set ──────────────────
list_t3 = Date[]
for yr in 2013:2025
    start_d = kdca_week_start(yr, 14); end_d = kdca_week_start(yr, 44) + Day(6)
    append!(list_t3, collect(start_d:Day(1):end_d))
end
valid_t3 = Date[]
for d in list_t3
    t2 = get(im_lookup, d, nothing)
    if !isnothing(t2) && year(d) == year(t2); push!(valid_t3, d); end
end
unique!(valid_t3)
possible_t4_set = Set{Date}()
for t3 in valid_t3, offset in 0:365
    t4_try = t3 - Day(offset)
    Date(2013,1,1) <= t4_try <= Date(2025,12,31) && push!(possible_t4_set, t4_try)
end
println("possible_t.4: $(length(possible_t4_set)) dates")

# ─── Parameters ──────────────────────────────────────────────────────────────
const t4_list = collect(Date(2013,1,2):Day(7):Date(2025,12,31))   # every Wednesday
println("t.4 list: $(length(t4_list)) dates")
const n_T_HM = 100; const n_T_MH = 30; const n_I_H = 500

# I_H (infectee incubation period): Gamma/LogNormal mixture, truncated to [1,500]
const d_IH_trunc = truncated(
    MixtureModel([Gamma(22.8197, 1.11405), LogNormal(5.78509, 0.140988)], [0.7423, 0.2577]),
    1.0, 500.0)

# ─── Core generation (single t4) ─────────────────────────────────────────────
function generate_SI_for_t4(t4::Date, alpha::Float64, beta::Float64, rng::AbstractRNG)
    t4 in possible_t4_set || return DataFrame()

    N_vec   = Float64[get_mos_N(t4 + Day(j)) for j in 0:365]
    cum_N   = cumsum(N_vec)
    Pr_T_HM = alpha .* N_vec .* exp.(-alpha .* cum_N)
    for j in 1:366
        Pr_T_HM[j] == 0.0 && continue
        t3_try = t4 + Day(j - 1)
        if t3_try < Date(2013,1,1) || t3_try >= Date(2026,1,1); Pr_T_HM[j] = 0.0; continue; end
        t2_try = get(im_lookup, t3_try, nothing)
        if isnothing(t2_try) || t2_try >= Date(2026,1,1); Pr_T_HM[j] = 0.0; continue; end
        year(t3_try) != year(t2_try) && (Pr_T_HM[j] = 0.0)
    end
    sum(Pr_T_HM) == 0.0 && return DataFrame()

    w_THM = Weights(Pr_T_HM)
    T_HM_samples = [sample(rng, 0:365, w_THM) for _ in 1:n_T_HM]
    t3_samples  = t4 .+ Day.(T_HM_samples)
    t2_samples  = [im_lookup[d]::Date for d in t3_samples]
    I_M_samples = [Dates.value(t2_samples[i] - t3_samples[i]) for i in 1:n_T_HM]

    n_4321  = n_T_HM * n_T_MH
    t4_4321 = fill(t4, n_4321)
    t3_4321 = repeat(t3_samples, inner=n_T_MH)
    t2_4321 = repeat(t2_samples, inner=n_T_MH)
    IM_4321 = repeat(I_M_samples, inner=n_T_MH)
    t1_4321 = Vector{Union{Date, Missing}}(missing, n_4321)

    for i in 1:n_T_HM
        u = t2_samples[i]; T_MH_max = 60 - I_M_samples[i]; base_idx = (i - 1) * n_T_MH
        T_MH_max < 0 && continue
        N_TMH   = Float64[get_mos_N(u + Day(j)) for j in 0:T_MH_max]
        Pr_T_MH = [(N_TMH[j+1] != 0.0 ? 1.0 : 0.0) * beta * exp(-beta * j) for j in 0:T_MH_max]
        sum(Pr_T_MH) == 0.0 && continue
        w_TMH = Weights(Pr_T_MH)
        for k in 1:n_T_MH
            t1_4321[base_idx + k] = u + Day(sample(rng, 0:T_MH_max, w_TMH))
        end
    end

    n_total  = n_4321 * n_I_H
    t4_final = repeat(t4_4321, inner=n_I_H)
    t3_final = repeat(t3_4321, inner=n_I_H)
    t2_final = repeat(t2_4321, inner=n_I_H)
    IM_final = repeat(IM_4321, inner=n_I_H)
    t1_final = repeat(t1_4321, inner=n_I_H)
    t0_final = Vector{Union{Date, Missing}}(missing, n_total)
    all_IH   = floor.(Int, rand(rng, d_IH_trunc, n_total))
    for row in 1:n_total
        ismissing(t1_final[row]) || (t0_final[row] = t1_final[row] + Day(all_IH[row]))
    end

    return DataFrame(t4 = t4_final, t3 = t3_final, t2 = t2_final,
                     t1 = t1_final, t0 = t0_final, I_M = IM_final)
end

# ─── Run 9 parameter combinations ────────────────────────────────────────────
for alpha_orig in ALPHA_LIST, beta_val in BETA_LIST
    alpha = alpha_orig / ALPHA_SCALE
    dn = joinpath(PROC, "SI_forward_a$(alpha_orig)_b$(beta_val)")
    mkpath(dn)
    if count(f -> endswith(f, ".jld2"), readdir(dn)) == length(t4_list)
        println("[α=$alpha_orig, β=$beta_val] SKIP (complete)"); continue
    end
    println("[α=$alpha_orig, β=$beta_val] Starting ($(length(t4_list)) dates)...")
    t_start = time(); save_lock = ReentrantLock()
    Threads.@threads for ind in 1:length(t4_list)
        t4 = t4_list[ind]; fn = joinpath(dn, "$(t4).jld2")
        isfile(fn) && continue
        dt = generate_SI_for_t4(t4, alpha, beta_val, MersenneTwister(ind))
        lock(save_lock) do
            jldopen(fn, "w"; iotype=IOStream, compress=true) do file; file["dt"] = dt; end
        end
    end
    println("[α=$alpha_orig, β=$beta_val] Done ($(round((time()-t_start)/60, digits=1)) min)")
end
println("\nAll 9 forward combinations completed.")
