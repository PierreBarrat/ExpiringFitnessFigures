### A Pluto.jl notebook ###
# v0.19.45

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
datdir = "data_trajectories_fixed_fitness.jl//"

# ╔═╡ e999242b-87b3-44a6-a264-9f1ffffc6101
files = open(datdir * "files.json", "r") do io
	JSON3.read(io, Dict)
end

# ╔═╡ c48c5983-1072-4793-a86b-a127bfbe71f3


# ╔═╡ 878b1b89-d74c-4fcc-be55-70d376f3b37e
md"# Probability of fixation"

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

# ╔═╡ 3cb788b1-ee3c-41c7-9663-72154a925aa3
begin
	svals = dat_all.s |> unique |> sort
	Δtvals = dat_all.Δt |> unique |> sort
	ρvals = dat_all.ρ |> unique |> sort

	Δt = 10
	fb = FrequencyBin(0.3, 0.05)
end

# ╔═╡ 7bb49e6c-5ec9-4333-bc94-8c8d867ced84
dat = begin
	dat_all[dat_all.Δt .== Δt, :];
end;

# ╔═╡ 5f687fbd-bbcd-4524-9eff-11f9493112c7
names(dat)

# ╔═╡ cffe4bb0-378d-4a0a-b51f-7fe9a3230958
always_below=true

# ╔═╡ d2aed60c-87dc-4f57-843a-605815a84aa4
function get_mean_trajectory(trajectories, fb; kwargs...)
	fmean, tmean, Zs = mean(trajectories, fb; always_below, kwargs...)
	return fmean, tmean
end

# ╔═╡ 7e5fcac2-844e-4ac0-97ca-706fc8092fcb
let 
	av_traj_df = select(
		dat, :ρ, :s, 
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
	
	local fmean, tmean, Zs = mean(trajectories, fb; always_below)
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
		title = "s = $(params.s), ρ = $(round(params.ρ, sigdigits=2)), s/ρ = $(round(params.s/params.ρ; sigdigits=2)), Δt = $(params.Δt)"
	)
end

# ╔═╡ 46587230-b730-4ff7-be1a-ab186e70e951
plt_dict = map(eachrow(dat)) do r
	(r.ρ, r.s) => plot_trajectories(r.trajectories, r, fb)
end |> Dict

# ╔═╡ 235b0068-51b6-458a-8cbd-4ffe82c35083
p1 = let
	plot(
		[plt_dict[ρ, s] for ρ in ρvals for s in svals]...,
		layout = grid(length(ρvals), length(svals)),
		size = (length(svals) * 600, length(ρvals) * 600),
		margin = 10mm,
		titlefontsize=18,
		legendfontsize=18,
		
	)
end

# ╔═╡ 642da2d5-d196-4fcf-936c-a9585a7df465
savefig(p1, datdir * "inertia_panel_f$(fb.f).png")

# ╔═╡ 54342047-a26e-4d2f-8f42-9f6bbb8f7816
function pfix_v_f(trajectories)
	fbs = [FrequencyBin(f, 0.05) for f in .2:.1:.8]
	pfix = 	map(fbs) do fb
		fixation_probability(T -> FT.inbin(T, fb; always_below), trajectories)
	end
	return (fbs, pfix)
end

# ╔═╡ 36004549-a98f-4663-a6db-18063d84a608
let
	pfix_df = select(
		dat, :ρ, :s, [:trajectories] => ByRow(pfix_v_f) => [:fb, :pfix]
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
		title = "s = $(params.s), ρ = $(round(params.ρ, sigdigits=2)), s/ρ = $(round(params.s/params.ρ; sigdigits=2)), Δt = $(params.Δt)"
	)
end

# ╔═╡ 8e892a2e-c51d-405c-aa3a-3036bcf36a2c
plt_dict_pfix = map(eachrow(dat)) do r
	(r.s, r.ρ) => plot_pfix(r.trajectories, r)
end |> Dict

# ╔═╡ c13f9cca-a534-4ecf-9b7a-0f8d4f825e55
p2 = let
	plot(
		[plt_dict_pfix[s, ρ] for s in svals for ρ in ρvals]...,
		layout = grid(length(svals), length(ρvals)),
		size = (length(ρvals) * 600, length(svals) * 600)
	)
end

# ╔═╡ 103ca2e1-cb08-406b-b29f-d5d53064f10f
plt_dict_pfix[0.03,0.03]

# ╔═╡ Cell order:
# ╠═eadd1d72-d1f5-11ed-098a-11e1b44a86a9
# ╠═a55b7fe0-798b-4eca-b795-ee1342721d08
# ╠═e999242b-87b3-44a6-a264-9f1ffffc6101
# ╠═c60dc463-ab42-4943-b6e7-c3c0b0efc3bd
# ╠═3cb788b1-ee3c-41c7-9663-72154a925aa3
# ╠═7bb49e6c-5ec9-4333-bc94-8c8d867ced84
# ╠═7e5fcac2-844e-4ac0-97ca-706fc8092fcb
# ╠═c48c5983-1072-4793-a86b-a127bfbe71f3
# ╠═46587230-b730-4ff7-be1a-ab186e70e951
# ╠═235b0068-51b6-458a-8cbd-4ffe82c35083
# ╠═642da2d5-d196-4fcf-936c-a9585a7df465
# ╟─878b1b89-d74c-4fcc-be55-70d376f3b37e
# ╠═36004549-a98f-4663-a6db-18063d84a608
# ╠═8e892a2e-c51d-405c-aa3a-3036bcf36a2c
# ╠═c13f9cca-a534-4ecf-9b7a-0f8d4f825e55
# ╠═103ca2e1-cb08-406b-b29f-d5d53064f10f
# ╟─f57bc80a-ab12-4911-bf98-3b3eefcb8ead
# ╠═6773d67d-2492-48b3-8c66-ee1199148606
# ╠═5f687fbd-bbcd-4524-9eff-11f9493112c7
# ╠═be559210-e29b-4369-b9f4-5232782b48df
# ╠═d2aed60c-87dc-4f57-843a-605815a84aa4
# ╠═ed3510c5-56f5-4404-9095-6d37cf6ad274
# ╠═64bf3d22-666d-49af-9ba0-70a3858bc5fb
# ╠═54342047-a26e-4d2f-8f42-9f6bbb8f7816
# ╠═cffe4bb0-378d-4a0a-b51f-7fe9a3230958
