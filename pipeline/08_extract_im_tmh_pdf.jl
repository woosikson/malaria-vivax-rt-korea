#!/usr/bin/env julia
#
# Step 8 — flatten the I_M + T_MH histogram to CSV (inputs for S1 Fig).
#   Reads the per-(α,β) normalised PDFs and means from step 05 and writes them
#   in long form.  day index 0..60 (day 0 = bin 1).
#
#   Input  : pipeline/processed/hist_IM_TMH.jld2  (from step 05)
#   Output : pipeline/processed/im_tmh_hist.csv   (alpha,beta,day,prob)
#            pipeline/processed/im_tmh_means.csv  (alpha,beta,mean)
#            → read by the R figure scripts from pipeline/processed/im_tmh_hist.csv, im_tmh_means.csv
#
# Run: julia pipeline/08_extract_im_tmh_pdf.jl

using JLD2, Printf
include(joinpath(@__DIR__, "_config.jl"))

src = joinpath(PROC, "hist_IM_TMH.jld2")
counts, means, alpha_list, beta_list = jldopen(src, "r") do f
    f["counts"], f["means"], f["alpha_list"], f["beta_list"]
end

open(joinpath(PROC, "im_tmh_hist.csv"), "w") do io
    println(io, "alpha,beta,day,prob")
    for a in alpha_list, b in beta_list
        p = counts["a$(a)_b$(b)"]
        for (d, prob) in enumerate(p)          # bin d=1 → day 0
            println(io, @sprintf("%s,%s,%d,%.10g", a, b, d - 1, prob))
        end
    end
end
open(joinpath(PROC, "im_tmh_means.csv"), "w") do io
    println(io, "alpha,beta,mean")
    for a in alpha_list, b in beta_list
        println(io, @sprintf("%s,%s,%.10g", a, b, means["a$(a)_b$(b)"]))
    end
end
println("Saved: im_tmh_hist.csv, im_tmh_means.csv")
