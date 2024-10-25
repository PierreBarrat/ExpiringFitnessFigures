### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 538982c4-615c-11ef-3f6e-e741e5c507a0
begin
	using Pkg; Pkg.activate("../../")
	using Chain
	using CSV
	using DataFrames
	using DataFramesMeta
	using DifferentialEquations
	using Distributions
	using FrequencyTrajectories
	using JSON3
	using LaTeXStrings
	using Measures
	using Parameters
	using Plots
	using PlutoUI
	using Random
	using StatsBase

	using DiscreteRandomWalk1D
	using FrequencyTrajectories
end

# ╔═╡ 223f967d-966c-43f1-b698-91483d1836c4
function pubfig(fnt_size=30; kwargs...)
    PLOTS_DEFAULTS = Dict(
        :markersize => 10,
        :linewidth => 3,
        :titlefontsize => fnt_size,
        :guidefontsize => fnt_size,
        :tickfontsize => fnt_size,
        :legendfontsize => fnt_size,
        :size => (1200,900),
        :gridlinewidth => 0.5,
        :gridalpha => 0.3,
        :framestyle => :box,
        :margin => 5mm,
        # :bottom_margin => 5mm,
        # :left_margin => 5mm,
    )

    for (k, x) in kwargs
        PLOTS_DEFAULTS[k] = x
    end

    return PLOTS_DEFAULTS
end

# ╔═╡ 0384a194-3d18-419c-9168-4e8c4ffd0a3a
begin
	plt_defaults = pubfig(24)
	Plots.default(; plt_defaults...)
end

# ╔═╡ a4bacac6-a461-4d0b-91d5-1267a1b5e0cd
md"## Expirng fitness trajectories"

# ╔═╡ d878c3c1-c256-4c52-aadc-a83f02a10294
readdir("../../")

# ╔═╡ 0a207907-ef3b-41b6-be11-ab75e11a6753
md"## flu trajecrories" 

# ╔═╡ b6cfe56d-fa84-48dd-9409-3666ffdaceb6
begin
	dat_h3n2 = open(io -> JSON3.read(io, Dict), "../h3n2_ha.json", "r")
	dat_sc2 = open(io -> JSON3.read(io, Dict), "../sc2.json", "r")
end

# ╔═╡ 848e3a83-3001-42a3-a080-6633ee6f2f63
md"## Utils"

# ╔═╡ 45ec10dc-895c-41c9-a7a2-00f1d99a54a8
function read_trajectories(file)
	return @chain file begin
		CSV.read(_, DataFrame) 
		FrequencyTrajectories.trajectories_from_dataframe
	end
end

# ╔═╡ 5134840b-624f-48c3-9d6a-c362b6fbd5f7
dat_exp_all = let
	datdir = joinpath(
		"predictability", 
		"trajectories_expiring_fitness",
		"data_trajectories_expiring_fitness.jl"
	)
	files = open(joinpath(datdir, "files.json"), "r") do io
		JSON3.read(io, Dict)
	end
	df = DataFrame(values(files))
	df.trajectories = map(df.trajectory_file) do f 
		read_trajectories(joinpath(datdir, f))
	end
	df
end;

# ╔═╡ 824b07d6-2a5a-43aa-94e8-f3d2e967c884
begin
	αvals = dat_exp_all.α |> unique |> sort
	ρvals = dat_exp_all.ρ |> unique |> sort
	svals = dat_exp_all.s |> unique |> sort
	Δtvals = dat_exp_all.Δt |> unique |> sort
end


# ╔═╡ 627c672f-a987-4af0-a854-05b547f7959c
ρvals

# ╔═╡ dc83b784-3899-46bb-89b5-6b9542f0f174
dat_exp = begin
	Δt = 10
	s = 0.03
	α = s
	subdat = @subset dat_exp_all begin
		:Δt .== Δt
		:α .== α
		:s .== s
		# :ρ .== ρ
	end
end

# ╔═╡ b9aaf2e6-a727-4f21-ad81-2df67a139a3c
function _plot_trajectories(trajectories, fb)
	tmin = -400
	tmax = 1200
	# fb = FT.FrequencyBin(.4, 0.03)
	
	local p = plot(
		ylim=(-0.02, 1.02), 
		xlim = (tmin, tmax),
		xlabel = "time",
		ylabel = "frequency",
		size = (900, 900),
		frame = :box
	)

	always_below = true
	# fmean, tmean, Zs = mean(trajectories, fb; K=2, always_below)
	for t in shuffle(trajectories)
		plot!(t, fb, label="", alpha = 0.8,)
	end
	
	# plot!(tmean, fmean, line=(:black, 3), label="")
	# plot!(collect(xlims(p)), [0.95, 0.95], label="", line=(:black, :dash))
	# plot!(collect(xlims(p)), [0.05, 0.05], label="", line=(:black, :dash))
	hline!([0.95], label="", line = (2, :black, :dash))
	hline!([0.05], label="", line = (2, :black, :dash))
	vline!([0.], line=(:grey, .5), label="")
	# plot!([tmin, 0], fb)
	
	return p
end

# ╔═╡ 259d177a-0164-452d-b505-61051da0ec57
function plot_trajectories(data, fb; N = Inf)
	# data is a DataFrame row
	trajectories = begin
		T = FrequencyTrajectories.filter(data.trajectories, fb; always_below=true)
		N > length(T) ? T : shuffle(T)[1:N]
	end
	p = _plot_trajectories(trajectories, fb)
	# vline!([params.Ne], line=(:black, :dash), label="Ne=$(params.Ne)")
	plot!(
		p,
		title = "ν = $(round(data.α, sigdigits=2)), ρ/s = $(round(data.ρ/data.s; sigdigits=2))"
	)
end

# ╔═╡ 69789b71-c4e5-425d-87f6-6e2afa60442b
clr_from_final(T, dx) = if T[end] > 1-dx
	:red
elseif dx < T[end] < 1-dx
	:black
else
	:blue
end

# ╔═╡ cf6f1904-dea0-43d1-b953-3a06a412d6f7
function plot_traj_flu(
	dat; 
	fbin = 0.5, sc2=false, write_ann=true, N = 10, kwargs...
)
	tvals = if !sc2
		dat["time_points"]
	else
		# length of trajectories
		# L = dat["trajectories"] |> values |> first |> values |> first |> length
		# @chain dat["time_points"] extrema range(_..., length=L) _.-50
		range(-7, 24) * 7
	end
	if !haskey(dat["trajectories"], string(fbin)) error() end

	p = plot(;
		xlabel = "days",
		ylabel = "frequency",
		title = "H3N2/HA",
		kwargs...
	)
	trajectories = @chain dat["trajectories"][string(fbin)] begin
		values
		collect
		shuffle
		getindex(1:N)
	end
	for T in trajectories
		Tc = convert(Vector{Union{Float64, Missing}}, T)
		# for (i, f) in enumerate(Tc)
		# 	if all(j -> Tc[j] < 0.05 || Tc[j] > .95, i:length(Tc))
		# 		Tc[i:end] .= missing
		# 		break
		# 	end
		# end
		# for (i, f) in Iterators.reverse(enumerate(Tc))
		# 	if all(j -> Tc[j] < 0.05 || Tc[j] > .95, 1:i)
		# 		Tc[1:i] .= missing
		# 		break
		# 	end
		# end
		color = clr_from_final(T, .2)
		plot!(tvals, Tc; label="", line=(.75, color))
	end

	# Tmean = mean(values(dat["trajectories"][string(fbin)]))
	# plot!(tvals, Tmean, line=(8, :black), label="")

	vline!([0.], line=(:grey, .5), label="")
	# hline!([0.95], label="", line = (2, :black, :dash))
	# hline!([0.05], label="", line = (2, :black, :dash))
	# plot!([tvals[1], 0], [fbin+.1, fbin+.1], line=(:black, .5), label="")
	# plot!([tvals[1], 0], [fbin, fbin]; line=(:black, .5), label="")

	(x_arr, x_ann) = sc2 ? (-50, -35) : (-320, -230)
	# plot!(
	# 	[x_arr, x_arr], [fbin, fbin+.1]; 
	# 	line=(3, :black), label="",  arrow = arrow(:both, .1, .2),
	# )
	# write_ann && annotate!(x_ann, fbin+0.05, text(L"x^\star", 36))
	
	# plot!([tvals[1], 0], [fbin-.05, fbin-0.05], line=(:black, .5), label="")
	p
end

# ╔═╡ ce7b9c6e-c797-42e8-a319-01e82a1556ca
let 
	Random.seed!(1)
	N = 10
	plt_expiring = map(eachrow(dat_exp)) do r 
		plot_trajectories(r, FrequencyBin(0.5, 0.025); N)
	end
	p = plot(
		plot_traj_flu(dat_h3n2; N), plt_expiring...;
		layout = grid(2, 2), size = (1600, 1600), dpi = 200,
	)
	savefig("../../figures/SI_figures/h3n2_vs_expiring_trajectories.png")
	p
end

# ╔═╡ 87ec6148-3abd-4fdc-a7a8-b0f076bceef2
plot_traj_flu(dat_h3n2; N=6)

# ╔═╡ Cell order:
# ╠═538982c4-615c-11ef-3f6e-e741e5c507a0
# ╠═223f967d-966c-43f1-b698-91483d1836c4
# ╠═0384a194-3d18-419c-9168-4e8c4ffd0a3a
# ╟─a4bacac6-a461-4d0b-91d5-1267a1b5e0cd
# ╠═5134840b-624f-48c3-9d6a-c362b6fbd5f7
# ╠═627c672f-a987-4af0-a854-05b547f7959c
# ╠═824b07d6-2a5a-43aa-94e8-f3d2e967c884
# ╠═dc83b784-3899-46bb-89b5-6b9542f0f174
# ╠═ce7b9c6e-c797-42e8-a319-01e82a1556ca
# ╠═d878c3c1-c256-4c52-aadc-a83f02a10294
# ╠═0a207907-ef3b-41b6-be11-ab75e11a6753
# ╠═b6cfe56d-fa84-48dd-9409-3666ffdaceb6
# ╠═87ec6148-3abd-4fdc-a7a8-b0f076bceef2
# ╠═848e3a83-3001-42a3-a080-6633ee6f2f63
# ╠═45ec10dc-895c-41c9-a7a2-00f1d99a54a8
# ╠═cf6f1904-dea0-43d1-b953-3a06a412d6f7
# ╠═b9aaf2e6-a727-4f21-ad81-2df67a139a3c
# ╠═259d177a-0164-452d-b505-61051da0ec57
# ╠═69789b71-c4e5-425d-87f6-6e2afa60442b
