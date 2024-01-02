### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 7e64bec1-9274-4ed3-8915-fcedd55db7b5
begin
	using Pkg; Pkg.activate("../")
	using Chain
	using DifferentialEquations
	using Distributions
	using Measures
	using Parameters
	using Plots
	using PlutoUI
	using Random
	using StatsBase

	using DiscreteRandomWalk1D
	using FrequencyTrajectories
	using PartialSweepSIR
end

# ╔═╡ a13a5083-5216-43ec-a469-fb4c28d0c7f0
begin
	# functions for individual plots
	include("./intermediate_notebooks/randomwalk_example_trajectories_asfunc.jl")
	include("./intermediate_notebooks/N_viruses_example_asfunc.jl")
	include("./intermediate_notebooks/ODE_expsel_asfunc.jl")
end

# ╔═╡ a15d04e6-338d-49dc-b867-1fb810d52d48
function pubfig(fnt_size=30; kwargs...)
    PLOTS_DEFAULTS = Dict(
        :markersize => 10,
        :linewidth => 5,
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

# ╔═╡ e3a280fe-889a-40de-aeb8-8a9e61636a18
begin
	plt_defaults = pubfig(20)
	Plots.default(; plt_defaults...)
end

# ╔═╡ c5eeae39-e578-43c9-93ea-6e25e12c3d6c
p_rw = let
	# good seeds I found for rw: 9, 15 
	seed_rw = 15
	p_rw = random_walk_example_plot(seed_rw)
	plot(
		p_rw; 
		ylim = (-0.025, 1.025),
		ann = (15, .92, ("D", 50)),
	)
end;

# ╔═╡ c3d7e947-5191-47a0-92d8-c6ad7a23ce6c
p_ode = let
	α = 0.05
	s0 = 0.05
	p = expiring_selection_plot(; α, s0)
	plot!(p.subplots[2], yticks = [0, s0])
	plot!(p.subplots[1], ylabel="")
	plot!(p, xlabel="", legendfontsize = 26, ann = (15, .85, ("C", 50)))
	plot!(p, legend = :topright)
end;

# ╔═╡ 70cc04b3-1d67-47f1-83d9-537c03fed930
expiring_selection_plot

# ╔═╡ 514cd78e-bc5d-4a31-ba96-66c7477098c7
p_freq, p_pop_stack = let
	# good seeds for long trajectories: 1, 5, 13, 21
	# good seeds for fixation: 3, 7, 10, 11, 12
	seed_pop_stack = 1
	ylabel = "Variant frequency"
	p_freq, p_pop_stack = pop_stack_plot(seed_pop_stack)
	(
		plot(
			p_freq; 
			title="", 
			xlabel="", 
			ylim = (-0.025, 1.025), 
			ylabel, 
			legendfontsize = 26,
			ann = (250, .92, ("A", 50))
		),
		plot(
			p_pop_stack; 
			title = "", ylim = (-0.025, 1.025), ylabel, ann = (250, .92, ("B", 50))
		)
	)
end

# ╔═╡ 8cf42471-7ed4-475b-a5d9-1b23b3912bdb
panel = let
	l = @layout [
		grid(2,2)
	]
	plot(
		p_freq, p_ode, p_pop_stack, p_rw; 
		layout=l, size = (2400, 1600), margin = 20mm, right_margin = 10mm, dpi=300,
	)
end

# ╔═╡ 0696369f-4bbe-4633-ba3d-8cdc2542acf7
savefig(panel, "../figures/expfit_randomwalk_panel.png")

# ╔═╡ 37aca45f-e245-4234-b2d7-b736716e87a6
savefig(joinpath(
	homedir(), 
	"Documents/BaleLabo/Notes/ExpiringFitness/figures/expfit_randomwalk_panel.png"
))

# ╔═╡ f94e0d7d-e3e8-4810-85ac-6f391488acb3
filler = plot(rand(10), title="FILLER");

# ╔═╡ Cell order:
# ╠═7e64bec1-9274-4ed3-8915-fcedd55db7b5
# ╠═a15d04e6-338d-49dc-b867-1fb810d52d48
# ╠═e3a280fe-889a-40de-aeb8-8a9e61636a18
# ╠═a13a5083-5216-43ec-a469-fb4c28d0c7f0
# ╠═c5eeae39-e578-43c9-93ea-6e25e12c3d6c
# ╠═c3d7e947-5191-47a0-92d8-c6ad7a23ce6c
# ╠═70cc04b3-1d67-47f1-83d9-537c03fed930
# ╠═514cd78e-bc5d-4a31-ba96-66c7477098c7
# ╠═8cf42471-7ed4-475b-a5d9-1b23b3912bdb
# ╠═0696369f-4bbe-4633-ba3d-8cdc2542acf7
# ╠═37aca45f-e245-4234-b2d7-b736716e87a6
# ╠═f94e0d7d-e3e8-4810-85ac-6f391488acb3
