#!/usr/bin/env julia
#
# Step 3 — weekly case reproduction number R_t (Binomial estimator) with 95% interval.
#   For each infector week t_j, each of its N_{t_j} cases is allocated to candidate
#   infectees t_i with probability p_ij = w_ij / Σ_k N_{t_k} w_kj (w = forward SI prob),
#   summed over Binomial(N_{t_i}, p_ij) draws, then divided by N_{t_j}.
#   10000 Monte-Carlo trials → 2.5/50/97.5% quantiles and mean.
#
#   Inputs : data/raw_inputs/malaria_kdca.csv
#            pipeline/processed/forward_infector_si_a{α}_b{β}.jld2  (from step 02)
#   Output : pipeline/processed/CI_R_case_a{α}_b{β}.csv  (date,alpha,beta,R_ll,R_med,R_mean,R_ul)
#            → read by the R figure scripts from pipeline/processed/CI_R_case_a{α}_b{β}.csv
#
#   NOTE: this step is stochastic (Binomial draws). A fixed seed is set so re-runs
#   are reproducible; values match the committed CSV to within Monte-Carlo error
#   (means stable to ~3 decimals over 10000 trials), which leaves the figures
#   visually identical.
#
# Run: julia --threads=auto pipeline/03_compute_R_case.jl

using CSV, DataFrames, Dates, Random, Distributions, Statistics, JLD2, Base.Threads
include(joinpath(@__DIR__, "_config.jl"))

println("Julia threads: $(Threads.nthreads())")
Random.seed!(20260616)

function load_case()::DataFrame
    df = DataFrame(CSV.File(CASE_CSV)); df.date = Date.(string.(df.date))
    return DataFrame(date = df.date, N = Int.(df.N))
end

function CI_R_case(df_case::DataFrame, alpha_label::Float64, beta::Float64)
    df_record_full = DataFrame(date_infectee = df_case.date, N = df_case.N)
    df_record = df_record_full[df_record_full.date_infectee .>= Date(2013,1,1) .&&
                               df_record_full.date_infectee .<  Date(2026,1,1), :]

    fn = joinpath(PROC, "forward_infector_si_a$(alpha_label)_b$(beta).jld2")
    si_dist_raw = jldopen(fn, "r"; iotype=IOStream) do file; file["si_dist"]; end

    n_dates = length(si_dist_raw)
    date_to_idx = Dict{Date, Int}()
    for i in 1:n_dates; date_to_idx[Date(si_dist_raw[i].date)] = i; end

    infector_dates = Date[]; infectee_dates = Union{Date, Missing}[]; si_weeks_vec = Int[]
    for i in 1:n_dates
        d_inf = Date(si_dist_raw[i].date); si_df = si_dist_raw[i].si_dist
        if nrow(si_df) == 0 || (nrow(si_df) == 1 && si_df.N[1] == 0)
            push!(infector_dates, d_inf); push!(infectee_dates, missing); push!(si_weeks_vec, 0)
        else
            for j in 1:nrow(si_df)
                sw = si_df.si_week[j]
                push!(infector_dates, d_inf); push!(si_weeks_vec, sw)
                push!(infectee_dates, d_inf + Day(7 * sw))
            end
        end
    end
    df_possible_si = DataFrame(date_infector = infector_dates, si_week = si_weeks_vec,
                               date_infectee = infectee_dates)

    get_si_dist(d::Date) = begin
        idx = get(date_to_idx, d, nothing)
        isnothing(idx) ? DataFrame() : si_dist_raw[idx].si_dist
    end
    date_to_N_full = Dict{Date, Int}(r.date_infectee => r.N for r in eachrow(df_record_full))

    infectee_to_infectors = Dict{Date, Vector{Date}}()
    for r in eachrow(df_possible_si)
        ismissing(r.date_infectee) && continue
        d_ee = r.date_infectee::Date
        haskey(infectee_to_infectors, d_ee) || (infectee_to_infectors[d_ee] = Date[])
        push!(infectee_to_infectors[d_ee], r.date_infector)
    end
    for (k, v) in infectee_to_infectors; infectee_to_infectors[k] = sort(unique(v)); end

    num_trials = 10000
    isod = df_record.date_infectee; n_isod = length(isod)
    R = zeros(num_trials, n_isod)
    println("  Computing R_t ($n_isod dates, $num_trials trials)..."); t_start = time()

    Threads.@threads for ind in 1:n_isod
        t_j = isod[ind]
        possible_t_i_all = df_possible_si[df_possible_si.date_infector .== t_j, :date_infectee]
        possible_t_i = unique(sort(skipmissing(possible_t_i_all) |> collect))
        filter!(d -> d < Date(2026,1,1), possible_t_i)
        N_tj = get(date_to_N_full, t_j, 0)
        if isempty(possible_t_i) || N_tj == 0; R[:, ind] .= 0.0; continue; end
        si_df_j = get_si_dist(t_j)
        for _ in 1:N_tj, t_i in possible_t_i
            si_week_val = div(Dates.value(t_i - t_j), 7)
            w_row = si_df_j[si_df_j.si_week .== si_week_val, :prob]
            isempty(w_row) && continue
            w = w_row[1]
            possible_t_k = get(infectee_to_infectors, t_i, Date[])
            denom = 0.0
            for t_k in possible_t_k
                N_tk = get(date_to_N_full, t_k, 0); si_df_k = get_si_dist(t_k)
                si_week_ki = div(Dates.value(t_i - t_k), 7)
                prob_row = si_df_k[si_df_k.si_week .== si_week_ki, :prob]
                isempty(prob_row) || (denom += N_tk * prob_row[1])
            end
            denom == 0.0 && continue
            N_ti = get(date_to_N_full, t_i, 0); N_ti == 0 && continue
            R[:, ind] .+= rand(Binomial(N_ti, w / denom), num_trials)
        end
        R[:, ind] ./= N_tj
    end
    println("  R_t done ($(round((time()-t_start)/60, digits=1)) min)")

    q(p) = [length(filter(!isnan, R[:, t])) == 0 ? 0.0 : quantile(filter(!isnan, R[:, t]), p) for t in 1:n_isod]
    R_mean = [length(filter(!isnan, R[:, t])) == 0 ? 0.0 : mean(filter(!isnan, R[:, t])) for t in 1:n_isod]
    return DataFrame(date = df_record.date_infectee,
        alpha = fill(alpha_label, n_isod), beta = fill(beta, n_isod),
        R_ll = Float64.(q(0.025)), R_med = Float64.(q(0.5)),
        R_mean = Float64.(R_mean), R_ul = Float64.(q(0.975)))
end

df_case = load_case()
println("N total = $(sum(df_case.N)) over $(nrow(df_case)) dates")
for alpha_label in ALPHA_LIST, beta in BETA_LIST
    out_fn = joinpath(PROC, "CI_R_case_a$(alpha_label)_b$(beta).csv")
    if isfile(out_fn); println("[α=$alpha_label, β=$beta] SKIP (output exists)"); continue; end
    println("[α=$alpha_label, β=$beta] Starting...")
    CSV.write(out_fn, CI_R_case(df_case, alpha_label, beta))
    println("[α=$alpha_label, β=$beta] Saved: $(basename(out_fn))"); flush(stdout)
end
println("\nAll 9 R_t computations completed.")
