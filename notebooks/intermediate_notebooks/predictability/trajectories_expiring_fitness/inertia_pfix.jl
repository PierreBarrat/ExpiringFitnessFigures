### A Pluto.jl notebook ###
# v0.19.37

using Markdown
using InteractiveUtils

# ╔═╡ eadd1d72-d1f5-11ed-098a-11e1b44a86a9
begin
	using Revise
	using Pkg; Pkg.activate("../../../../")
	using Chain
	using CSV
	using DataFrames
	using FrequencyTrajectories
	using JSON3
	using Measures
	using Plots
	using Random
	using StatsBase
end

# ╔═╡ a55b7fe0-798b-4eca-b795-ee1342721d08
datdir = "data_trajectories_expiring_fitness.jl/"

# ╔═╡ e999242b-87b3-44a6-a264-9f1ffffc6101
files = open(datdir * "files.json", "r") do io
	JSON3.read(io, Dict)
end

# ╔═╡ c06c8d03-1691-4962-a054-6c168646b436
md"""
## Writing statistics to files
"""

# ╔═╡ 46587230-b730-4ff7-be1a-ab186e70e951
# ╠═╡ disabled = true
#=╠═╡
plt_dict = map(eachrow(dat)) do r
	(r.α, r.ρ) => plot_trajectories(r.trajectories, r, fb)
end |> Dict
  ╠═╡ =#

# ╔═╡ 235b0068-51b6-458a-8cbd-4ffe82c35083
#=╠═╡
p1 = let
	plot(
		[plt_dict[α, ρ] for α in αvals for ρ in ρvals]...,
		layout = grid(length(αvals), length(ρvals)),
		size = (length(ρvals) * 600, length(αvals) * 600),
		margin = 10mm,
		titlefontsize=18,
		legendfontsize=18,
		
	)
end
  ╠═╡ =#

# ╔═╡ 642da2d5-d196-4fcf-936c-a9585a7df465
#=╠═╡
savefig(p1, datdir * "inertia_panel_f$(fb.f).png")
  ╠═╡ =#

# ╔═╡ 60dc2d79-ca84-41f9-aab4-8388cec6f1cb
md"# Probability of fixation"

# ╔═╡ 063cb3a3-1a08-4b3c-b407-57541a93581c


# ╔═╡ f57bc80a-ab12-4911-bf98-3b3eefcb8ead
md"## Functions"

# ╔═╡ 6773d67d-2492-48b3-8c66-ee1199148606
function read_trajectories(file)
	return @chain begin
		file
		datdir * _
		CSV.read(_, DataFrame) 
		FrequencyTrajectories.trajectories_from_dataframe
	end
end

# ╔═╡ c60dc463-ab42-4943-b6e7-c3c0b0efc3bd
dat_all = let
	df = DataFrame(values(files))
	df.trajectories = map(f -> read_trajectories(f), df.trajectory_file)
	df
end;

# ╔═╡ 9757d8d3-d318-4b38-a322-ed1261f351f2
begin
	αvals = dat_all.α |> unique |> sort
	ρvals = dat_all.ρ |> unique |> sort
	svals = dat_all.s |> unique |> sort
	Δtvals = dat_all.Δt |> unique |> sort
end

# ╔═╡ 5804c6d3-72c0-4435-8a55-47db95092866
begin
	@assert length(svals) == 1 
	s0 = svals[1]
	Δt = 10
	fb = FrequencyBin(0.5, 0.01)
end

# ╔═╡ 7bb49e6c-5ec9-4333-bc94-8c8d867ced84
dat = begin
	dat_all[dat_all.Δt .== Δt, :];
end;

# ╔═╡ 549a0440-5e19-4359-890b-2879c1c1bdbb
names(dat_all)

# ╔═╡ 5f687fbd-bbcd-4524-9eff-11f9493112c7
names(dat)

# ╔═╡ 14e7e798-adcf-4d3e-bbab-0d09eb73874a


# ╔═╡ cffe4bb0-378d-4a0a-b51f-7fe9a3230958
always_below=true

# ╔═╡ 3f5f64e2-ed11-47a7-a7b1-a4ceabde628e
function get_mean_trajectory(trajectories, fb; kwargs...)
	fmean, tmean, Zs = mean(trajectories, fb; always_below, kwargs...)
	return fmean, tmean
end

# ╔═╡ 87c662b4-1eb5-4247-9632-e3a9297de5dc
let 
	av_traj_df = select(
		dat, :α, :ρ, :s, 
		[:trajectories] => ByRow(t -> get_mean_trajectory(t, fb; K=1/2)) => [:mf, :mt]
	)
	CSV.write("data_avtraj.csv", av_traj_df)
end

# ╔═╡ ed3510c5-56f5-4404-9095-6d37cf6ad274
function _plot_trajectories(trajectories, fb)
	tmin = -400
	tmax = 1200
	# fb = FT.FrequencyBin(.4, 0.03)
	
	local p = plot(
		ylim=(-0.02, 1.02), 
		xlim = (tmin, tmax),
		xlabel = "Time",
		ylabel = "Frequency",
		size = (900, 900),
		frame = :box
	)
	
	local fmean, tmean, Zs = mean(trajectories, fb; K=2, always_below)
	for t in shuffle(FrequencyTrajectories.filter(trajectories, fb; always_below))
		plot!(t, fb, label="", alpha = 0.25,)
	end
	
	plot!(tmean, fmean, line=(:black, 3), label="")
	plot!(collect(xlims(p)), [0.95, 0.95], label="", line=(:black, :dash))
	plot!(collect(xlims(p)), [0.05, 0.05], label="", line=(:black, :dash))
	plot!([tmin, 0], fb)
	
	return p
end

# ╔═╡ be559210-e29b-4369-b9f4-5232782b48df
function plot_trajectories(trajectories, params, fb)
	p = _plot_trajectories(trajectories, fb)
	# vline!([params.Ne], line=(:black, :dash), label="Ne=$(params.Ne)")
	plot!(
		p,
		title = "s = $(params.s), ρ = $(round(params.ρ, sigdigits=2)), α = $(round(params.α, sigdigits=2)), s/ρ = $(round(params.s/params.ρ; sigdigits=2))"
	)
end

# ╔═╡ 54342047-a26e-4d2f-8f42-9f6bbb8f7816
function pfix_v_f(trajectories)
	fbs = [FrequencyBin(f, 0.01) for f in .2:.1:.8]
	pfix = 	map(fbs) do fb
		fixation_probability(T -> FT.inbin(T, fb; always_below), trajectories)
	end
	return (fbs, pfix)
end

# ╔═╡ 71d8a61b-dd96-4119-965e-6e0b8cbdee45
let
	pfix_df = select(
		dat, :α, :ρ, :s, :mβ, :β2, [:trajectories] => ByRow(pfix_v_f) => [:fb, :pfix]
	) 
	select!(pfix_df, Not(:fb), :fb => ByRow(col -> map(fb -> fb.f, col)) => :f)
	CSV.write("data_pfix.csv", pfix_df)
end;

# ╔═╡ 06b5f62e-a53a-4a38-a153-a806648bb242
let
	pfix_df = select(
		dat, :α, :ρ, :s, [:trajectories] => ByRow(pfix_v_f) => [:fb, :pfix]
	) 
	select!(pfix_df, Not(:fb), :fb => ByRow(col -> map(fb -> fb.f, col)) => :f)
	CSV.write("data_pfix.csv", pfix_df)
end

# ╔═╡ 64bf3d22-666d-49af-9ba0-70a3858bc5fb
function plot_pfix(trajectories, params)
	fbs, pfix = pfix_v_f(trajectories)

	p = plot([fb.f for fb in fbs], pfix, label="", marker=:o)
	plot!([0,1], [0,1], line=(:black, :dash), label="")
	plot!(
		p,
		title = "s = $(params.s), ρ = $(round(params.ρ, sigdigits=2)), α = $(round(params.α; sigdigits=2))",
	)
end

# ╔═╡ 8e892a2e-c51d-405c-aa3a-3036bcf36a2c
plt_dict_pfix = map(eachrow(dat)) do r
	(r.α, r.ρ) => plot_pfix(r.trajectories, r)
end |> Dict

# ╔═╡ c13f9cca-a534-4ecf-9b7a-0f8d4f825e55
p2 = let
	plot(
		[plt_dict_pfix[a, ρ] for a in αvals for ρ in ρvals]...,
		layout = grid(length(αvals), length(ρvals)),
		size = (length(ρvals) * 400, length(αvals) * 400)
	)
end

# ╔═╡ 4109d6ab-f275-4d1b-b159-a745da80360e
length(plt_dict_pfix)

# ╔═╡ Cell order:
# ╠═eadd1d72-d1f5-11ed-098a-11e1b44a86a9
# ╠═a55b7fe0-798b-4eca-b795-ee1342721d08
# ╠═e999242b-87b3-44a6-a264-9f1ffffc6101
# ╠═c60dc463-ab42-4943-b6e7-c3c0b0efc3bd
# ╠═9757d8d3-d318-4b38-a322-ed1261f351f2
# ╠═5804c6d3-72c0-4435-8a55-47db95092866
# ╠═7bb49e6c-5ec9-4333-bc94-8c8d867ced84
# ╟─c06c8d03-1691-4962-a054-6c168646b436
# ╠═71d8a61b-dd96-4119-965e-6e0b8cbdee45
# ╠═87c662b4-1eb5-4247-9632-e3a9297de5dc
# ╠═46587230-b730-4ff7-be1a-ab186e70e951
# ╠═235b0068-51b6-458a-8cbd-4ffe82c35083
# ╠═642da2d5-d196-4fcf-936c-a9585a7df465
# ╟─60dc2d79-ca84-41f9-aab4-8388cec6f1cb
# ╠═06b5f62e-a53a-4a38-a153-a806648bb242
# ╠═549a0440-5e19-4359-890b-2879c1c1bdbb
# ╠═8e892a2e-c51d-405c-aa3a-3036bcf36a2c
# ╠═063cb3a3-1a08-4b3c-b407-57541a93581c
# ╠═c13f9cca-a534-4ecf-9b7a-0f8d4f825e55
# ╠═4109d6ab-f275-4d1b-b159-a745da80360e
# ╟─f57bc80a-ab12-4911-bf98-3b3eefcb8ead
# ╠═6773d67d-2492-48b3-8c66-ee1199148606
# ╠═5f687fbd-bbcd-4524-9eff-11f9493112c7
# ╠═be559210-e29b-4369-b9f4-5232782b48df
# ╠═3f5f64e2-ed11-47a7-a7b1-a4ceabde628e
# ╠═ed3510c5-56f5-4404-9095-6d37cf6ad274
# ╠═64bf3d22-666d-49af-9ba0-70a3858bc5fb
# ╠═54342047-a26e-4d2f-8f42-9f6bbb8f7816
# ╠═14e7e798-adcf-4d3e-bbab-0d09eb73874a
# ╠═cffe4bb0-378d-4a0a-b51f-7fe9a3230958
