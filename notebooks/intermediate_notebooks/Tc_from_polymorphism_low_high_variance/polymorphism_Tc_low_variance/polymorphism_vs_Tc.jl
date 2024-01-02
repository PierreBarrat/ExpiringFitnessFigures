### A Pluto.jl notebook ###
# v0.19.30

using Markdown
using InteractiveUtils

# ╔═╡ 5127cf4e-eb40-11ed-0a15-53a1acf43aef
begin
	using Revise
	using Pkg; Pkg.activate("../../../../")
	using Chain
	using CSV
	using DataFrames
	using DelimitedFiles
	using Distributions
	using FrequencyTrajectories
	using JSON3
	using Measures
	using Plots
	using QuadGK
	using Random
	using StatsBase
	using SpecialFunctions
end

# ╔═╡ 5ffd9c4e-4807-43aa-943a-4c8c3f9b022a
datdir = "data_trajectories_random_beta.jl/"

# ╔═╡ 2713455c-df95-49c7-925a-9d92b91719a5
files = open(datdir * "files.json", "r") do io
	JSON3.read(io, Dict)
end

# ╔═╡ 0e6415c0-b44d-4c8d-9b76-806d9f43926f


# ╔═╡ 9ac16620-e45c-4218-9d7f-c7cca89e66e1
md"# Functions"

# ╔═╡ 9ff2d96e-4637-4699-a4a2-47f4d9887fa5
"""
	get_β_distribution(mβ, β2)

Return a Beta distribution with mean `mβ` and second moment `β2`.
"""
function get_β_distribution(mβ, β2)
	b = (mβ - β2)/(β2 - mβ^2)*(1-mβ)
	a = mβ/(1-mβ)*b
	return Beta(a,b)
end

# ╔═╡ c9a65c7c-2b3b-403a-9141-6bc9fb316512
function fraction_covered_by_sweeps(mβ, β2, ρ, α)
	B = get_β_distribution(mβ, β2)
	return - 3 * ρ / α / mean(β -> log(1-β), rand(B, 1000))
end

# ╔═╡ fae29d40-0bd5-43cb-ae2d-e264df34719e
function frac_overlap(mβ, β2, ρ, α; K = 2)
	# parameters of the B distribution
	b = (mβ - β2)/(β2 - mβ^2)*(1-mβ)
	a = mβ/(1-mβ)*b
	# probability that s < ρ for a given sweep
	integrand(t) = ρ*exp(-ρ*t)*beta_inc(a, b, 1-exp(-K/α/t))[1]
	return quadgk(t -> integrand(t), 0, 5/ρ, rtol=1e-5)[1]
end

# ╔═╡ 10987b07-70e0-4a19-9173-83cf487c08ff
dat = let
	df = DataFrame(values(files))
	df.poly = map(df.diversity_file) do f
		X = CSV.read(datdir * f, DataFrame)
		X.polymorphism / 2
	end
	df.tvals =  map(df.diversity_file) do f
		X = CSV.read(datdir * f, DataFrame)
		X.t
	end

	df.sweep_fraction = map(eachrow(df)) do r
		fraction_covered_by_sweeps(r.mβ, r.β2, r.ρ, r.α)
	end
	df.frac_overlap = map(eachrow(df)) do r
		frac_overlap(r.mβ, r.β2, r.ρ, r.α)
	end
	# this proxy does not work very well
	# df.β2_eff = map(eachrow(df)) do r
	# 	β2_effective(r.mβ, r.β2, r.ρ, r.α)
	# end
	
	df.Tc_theory = map(r -> 1/r.ρ/r.β2, eachrow(df))
	df.Tc_measured = map(eachrow(df)) do r
		x = mean(r.poly[1:end])
		x/r.μ_neutral
	end
	
	sort(df, [:ρ, :β2])
end;

# ╔═╡ 22a5f98d-6b07-4d8c-bafa-dfeed71e9749
begin
	# export data used for panel
	local cols = vcat(
		map(x -> x => x, [:mβ, :β2, :ρ, :α, :frac_overlap]), 
		[:Tc_theory => :Ne, :Tc_measured => :av_poly],
	)
	df_export = select(dat, cols)
	CSV.write("data_panel.csv", df_export)
end

# ╔═╡ 94813cd3-210a-4b77-aac7-f68e1dc3a1ba
begin
	# Constants
	α = dat.α[1]
	μ_neutral = dat.μ_neutral[1]
	L_sel = dat.L[1]
	L_neutral = dat.L_neutral[1]
	L = L_sel + L_neutral
	neutral_sites = (1+2*L_sel):2*L
end

# ╔═╡ 9b6ba4bb-a5e6-4f19-9da3-bb9690c6c0aa
parameters = map(r -> (mβ=r.mβ, β2=r.β2, ρ=r.ρ), eachrow(dat)) |> sort |> unique

# ╔═╡ 6a5d6c8f-618d-4cb4-a519-31042c459f5f
begin
	mβ_vals = dat.mβ |> unique |> sort
	β2_vals = dat.β2 |> unique |> sort
	ρ_vals = dat.ρ |> unique |> sort
end

# ╔═╡ 87b9ece6-eb02-407b-b7ee-161e11932c32
dist_plts = map(zip(mβ_vals, β2_vals)) do (mβ, β2)
	B = get_β_distribution(mβ, β2)
	xvals = range(0, 1, length=1000)
	p = plot(xlabel="β", frame=:box, title="β distribution")
	plot!(xvals, map(x -> pdf(B, x), xvals), label="", line=(:black))
	p
end

# ╔═╡ 51cdd66e-99d3-42f2-8dba-941b54ad6e34
function Pτ(β, B, x0, ρ, α)
	τ = -1/α * (log(β) - log(x0)) / log(1-β)
	# τ = -1/α / log(1-β)
	exp(-ρ *τ)
end

# ╔═╡ 72eb54ac-7424-496f-abe2-f43c687e8e6b
function β2_effective(mβ, β2, ρ, α)
	# parameters of the B distribution
	b = (mβ - β2)/(β2 - mβ^2)*(1-mβ)
	a = mβ/(1-mβ)*b
	B = Beta(a, b)
	# probability that s < ρ for a given sweep
	x0 = 0.01
	I = quadgk(2*x0, 1-2*x0, rtol=1e-5) do β
		pdf(B, β) * β^2 * Pτ(β, B, x0, ρ, α)
	end

	return I[1]
end

# ╔═╡ 83eb61b1-622b-4b1d-a1bf-3a7a0a750cb5
function sweep_time(mβ, β2, α)
	B = get_β_distribution(mβ, β2)
	return @chain begin
		rand(B, 10_000) # β sample
		map(β -> log(1-β), _)
		filter(isfinite, _)
		mean
		- 3 / α / _
	end
end

# ╔═╡ 04a1b0ee-233f-4c94-a150-c95ffd10c868
function plot_polymorphism(r)
	idx = 1:10:length(r.poly)
	plot(r.tvals[idx], r.poly[idx] / r.μ_neutral, label="", alpha=0.3)
	hline!([r.Tc_theory], label="Ne", line=(:red, 3, 0.9))

	hline!([r.Tc_measured], label="<x(1-x)>/μ", color=1, line=(:dashdot, 4))
	# hline!([1/r.ρ/r.β2_eff])
	st = round(sweep_time(r.mβ, r.β2, r.α), sigdigits=2)
	dt = round(Int, 1/r.ρ)
	fOVL = round(r.frac_overlap, sigdigits=2)

	xticks = map(range(extrema(r.tvals[idx])..., length=3)) do x
		round(x, sigdigits=2) 
	end
	
	plot!(
		xlabel= r.mβ == maximum(mβ_vals) ? "time" : "",
		ylabel = r.ρ == minimum(ρ_vals) ? "x(1-x)/μ" : "",
		# title = "1/ρ = $(dt), sweep time ~ $(st)",
		title = """ρ=$(round(r.ρ, sigdigits=2)), <β> = $(r.mβ)
		sweep overlap ~ $(fOVL)
		""",
		# ylim = (-100, 6*Ne),
		tickfontsize = 6,
		frame=:box,
		xticks=xticks,
		# formatter = x -> format(x, precision=2)
	)
end

# ╔═╡ bd8b7c82-c16f-4fbf-ac92-6e213580a641
plot_polymorphism(dat[4,:])

# ╔═╡ 21c7431a-3b6d-41f0-b5fb-e989bad2259a
plts = map(Iterators.product(ρ_vals, β2_vals)) do (ρ, β2)
	df = dat[dat.β2 .== β2 .&& dat.ρ .== ρ, :]
	plot_polymorphism(first(eachrow(df)))
end;

# ╔═╡ 88229ffb-bdb1-4691-9f1a-2e8eccec2954
let
	p = plot(
		# plts...;
		vcat(plts, reshape(dist_plts, 1, 3))...;
		layout = (length(β2_vals), length(ρ_vals)+1),
		size = (500*length(ρ_vals), 500*(length(β2_vals)+1)),
		bottommargin = 10mm,
		topmargin = 10mm,
		leftmargin = 20mm,
		guidefontsize = 16,
		legendfontsize = 16,
		tickfontsize = 16, 
		
	)
	savefig("panel_diversity_v_time.png")
	p
end

# ╔═╡ Cell order:
# ╠═5127cf4e-eb40-11ed-0a15-53a1acf43aef
# ╠═5ffd9c4e-4807-43aa-943a-4c8c3f9b022a
# ╠═2713455c-df95-49c7-925a-9d92b91719a5
# ╠═10987b07-70e0-4a19-9173-83cf487c08ff
# ╠═22a5f98d-6b07-4d8c-bafa-dfeed71e9749
# ╠═94813cd3-210a-4b77-aac7-f68e1dc3a1ba
# ╠═04a1b0ee-233f-4c94-a150-c95ffd10c868
# ╠═bd8b7c82-c16f-4fbf-ac92-6e213580a641
# ╠═9b6ba4bb-a5e6-4f19-9da3-bb9690c6c0aa
# ╠═6a5d6c8f-618d-4cb4-a519-31042c459f5f
# ╠═21c7431a-3b6d-41f0-b5fb-e989bad2259a
# ╠═87b9ece6-eb02-407b-b7ee-161e11932c32
# ╠═88229ffb-bdb1-4691-9f1a-2e8eccec2954
# ╠═0e6415c0-b44d-4c8d-9b76-806d9f43926f
# ╠═9ac16620-e45c-4218-9d7f-c7cca89e66e1
# ╠═9ff2d96e-4637-4699-a4a2-47f4d9887fa5
# ╠═c9a65c7c-2b3b-403a-9141-6bc9fb316512
# ╠═fae29d40-0bd5-43cb-ae2d-e264df34719e
# ╠═72eb54ac-7424-496f-abe2-f43c687e8e6b
# ╠═51cdd66e-99d3-42f2-8dba-941b54ad6e34
# ╠═83eb61b1-622b-4b1d-a1bf-3a7a0a750cb5
