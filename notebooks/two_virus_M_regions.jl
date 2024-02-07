### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 4982e04e-6501-11ed-132e-ad3e395233da
begin
	using Pkg; Pkg.activate("../")
	using Chain
	using Measures
	using Parameters
	using Plots
	using PlutoUI
	using PartialSweepSIR
	using StatsBase
end

# ╔═╡ c2907ee7-3c4e-4b6b-946b-85971c5abbff
Plots.backend()

# ╔═╡ 6ab4e0ec-8ff1-4415-bcb8-5e5fbaf45bfc
md"""
Use like this:
```
font_size = 20
plt_defaults = pubfig(font_size) # pubfig() defaults to 30
Plots.default(; plt_defaults...)
```
"""

# ╔═╡ 5e0b0485-e3dc-401a-818f-bcd8924995b5
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

# ╔═╡ f905453d-2b68-40fd-9950-a7165933e73f
begin
	plt_defaults = pubfig(24)
	Plots.default(; plt_defaults...)
end

# ╔═╡ db663438-8523-47ac-a783-88ff59a7ab1d
md"# Setup"

# ╔═╡ c6c407d9-fea7-4e79-97ff-6acb88f4356d
md"""
Two sliders: 
- `M` controls the number of regions;
- `c` controls the connectivity between regions, with $0 \leq c \leq 1/M$

The `C` matrix can be accesssed in `params.C`
"""

# ╔═╡ 9c929838-c809-42c2-81b3-546a9de52e94
Ms = let
	_Ms = @bind M Slider([1,2,3,4,5,10,20], show_value=true, default=10)
	md"M = $(_Ms)"
end

# ╔═╡ 943ee06f-d4e2-40fe-9764-6c0a0b7ef0b0
I0 = 1e-4 # initial level of the variant

# ╔═╡ a7fc518a-f6d3-4e93-bac5-02d946c768f8
md"Cross immunity matrix for the special region $i=1$."

# ╔═╡ c4da46f5-8178-4ae5-b865-db025a365da0
K_special = begin
	b_special = .8
	f_special = .65
	[1 b_special; f_special 1]
end

# ╔═╡ 732be632-4b69-4434-bf51-33a219632fe2
β = (1-f_special)/(2-b_special-f_special)

# ╔═╡ 81cbea6e-9a97-4e7e-9528-f0830c7ae58f
md"Cross immunity matrix for other regions $i>1$."

# ╔═╡ 78b03821-37d4-4a9d-8ff0-349e12f66686
K0 = begin
	x = .99
	b = x
	f = x
	[1 b; f 1]
end

# ╔═╡ b7dd719a-118e-4809-84da-609735534835
md"""
Setting up groups in a vector `R`: 
- `R[1]` has a special cross-immunity matrix
- other groups have `K0`
"""

# ╔═╡ a62f2142-7d76-4a15-9566-54115aa98e08
one_group, many_groups = let
	S0 = .4
	I0 = 1e-6
	C0 = 0
	v1 = PSS.Virus(; S=S0, I=I0, C=C0) # wt virus
	v2 = PSS.Virus(; S=S0, I=0, C=0.) # mutant: initially not here

	special_group = PSS.Region(; viruses=[v1,v2], K=K_special)
	R = [PSS.Region(; viruses=[v1,v2], K=K0) for m in 1:M]
	R[1] = special_group
	
	[special_group], R
end

# ╔═╡ 51c77696-d4ae-466e-a19a-cf9278c7d3a5
md"# Simulation"

# ╔═╡ 559e0978-3bb1-407f-b27b-4e31dc9df04e
md"Simulating for time `T` to get the initial state, and then introducing the mutant in all regions."

# ╔═╡ 64c6f22b-9717-46c4-93cf-e1eafd107931
md"## Quick note on indexing"

# ╔═╡ 1d3084e3-a18e-4850-b424-a7a32f4780b3
md"""
The type `PartialSweepSIR.SIRState`, like the variables `state_init` and `state_final`, can be indexed in the following way: 
- `state[i, a, X]` where `X` is in `[:S, :I, :C, :R]` returns the size of compartment `X` in region `i` for virus `a`. `i` and `a` can also be ranges, *e.g.* `state[1, 1:2, :I]`.
- `state[:, a, X]` returns a vector of size `M` with the size of compartment `X` in all regions $i\in[1,M]$
- `state[i, :, X]`: same, but for all viruses
"""

# ╔═╡ 7aa5f567-8af0-4b40-8669-f959da124826
md"""
The type `PartialSweepSIR.Solution`, like the variable `sol`, can be used like this: 
- `sol(t)` returns the state at time `t` in the form of an `SIRState` object. The paragraph above shows how to index it.
- `sol[tvals, i, a, X]` returns the vector of the size of compartment `X` in region `i` for virus `a` and for all values in `tvals`.  
"""

# ╔═╡ ffebafd7-5abf-461f-af8d-030448f77f44
md"# Figures"

# ╔═╡ 29fa649f-7400-4099-8163-75a3d12362c9
backend()

# ╔═╡ 6c3a3b73-34b2-4643-9170-ed644e363725
md"## Panel"

# ╔═╡ d65daa41-6477-4eef-8193-364ee4cc21f6
md"## Frequency"

# ╔═╡ d1910dd1-3e5e-4d1c-92cb-529b9a00c356
default(:fontfamily)

# ╔═╡ 8f365998-cb7e-4000-a5c5-516a2b1ac814
function logistic_growth(tvals, s, x0)
	return map(tvals) do t
		exp(s*t) / (1/x0 - 1 + exp(s*t))
	end
end

# ╔═╡ 38e8868e-43ff-4f8b-88fb-07930638c0dc
md"## Infectious"

# ╔═╡ d4d20d93-d54b-4fb8-a334-dde0a6ef5684
md"## Susceptibles"

# ╔═╡ cad2211c-3f9f-43e7-8098-17c07dc4c12f
S_yrange = (0.25, 0.44)

# ╔═╡ 4640dbe6-89b5-4590-ac21-914226ca2097
md"# Misc"

# ╔═╡ c79fe879-a32a-4391-8df1-5b9d9eeb75f5
function log_range(low, high, length)
	@assert low > 0 && high > low
	exp.(range(log(low), log(high); length))
end

# ╔═╡ 98f8262c-28e4-4be1-91e9-6ac96bbb72a8
Cs = let
	c_values = vcat([0], log_range(1/M*1e-3, 1/M, 10))
	_Cs = @bind c Slider(c_values, show_value=true, default = 1/M/10)
	md"c = $(_Cs)"
end

# ╔═╡ f431ba40-9742-4898-92f0-33c3b8e12a1c
params_one, params_many, params_many_varC = let
	N = 2
	α = 3
	γ = .005
	(
		PSS.Parameters(; N, M=1, α, γ, c=1), 
		PSS.Parameters(; N, M, α, γ, c=1/M),
		PSS.Parameters(; N, M, α, γ, c),
	)
end;

# ╔═╡ e38f872f-f9ec-45cf-953c-42a18aa0bab3
@unpack α, γ, δ = params_one;

# ╔═╡ d2b8e952-37e1-4332-a321-706bd589a62d
begin
	Seq_special = δ/(δ + f_special*(α-δ))
	Seq_other = δ/(δ + x*(α-δ))

	growth_rate_one = α*Seq_special - δ
	growth_rate_many = α*((1-1/M)*Seq_other + Seq_special/M) - δ
end

# ╔═╡ c2a159f4-cd25-4525-a270-2b3065a0dc0e
T = 5/params_one.γ # simulation time

# ╔═╡ fd0d5424-074d-4fee-a41c-a80e61591348
state_init_many, sol_init_many, sol_many = let
	state = PSS.SIRState(; regions=many_groups, parameters=params_many)
	# simulating for  time T without the variant to get initial state
	sol_init = PSS.simulate(state, (0, T))
	# introduce the variant
	state_init = PSS.set_infected(sol_init(T), 2, I0)
	# re-simulate for time T
	sol = PSS.simulate(state_init, (0, T))

	state_init, sol_init, sol
end;

# ╔═╡ 0e026c1f-a6c9-49d5-8cf2-23d0b33a217f
state_init_many.regions[2].viruses[2]

# ╔═╡ d33785f7-96b6-4c26-92a7-82314b8114e1
p_freq_many = let
	# freq. plot
	tvals = range(0, T, length=100)
	f_R1 = PSS.frequency(sol_many, tvals, 1, 2)
	# f_R2 = PSS.frequency(sol, tvals, 2, 2)
	f_R = PSS.frequency(sol_many, tvals, 2)

	p = plot(
		legend = :bottomright,
		xlabel = "",
		ylabel = "",
		ylim = (-0.01, 1.01),
	)

	plot!(tvals, f_R, line=(:blue), label="")
		
	hline!([β]; label="", line=(:gray, 10, .5))

	annotate!(T*.95, β+.1, ("β", 40))
	annotate!(T*.15, .17, ("M=$M immune groups", 30, :left,  "Helvetica Bold"))
	annotate!(T*.15, .07, ("Fast mixing", 30, :left, "Helvetica Bold"))
	# plot logistic growth for comparison
	plot!(tvals, logistic_growth(
		tvals, growth_rate_many, f_R[1]);
		label = "", line=(:black, :dash, 2),
	)
	
	p
end

# ╔═╡ 2f8e4f1f-2da0-4bfa-b63c-069082180b6d
pI_many = let
	g = :I
	tvals = range(0, T, length=400)
	
	X_wt = map(t -> mean(sol_many[t, 1:M, 1, g]), tvals)
	X_m = map(t -> mean(sol_many[t, 1:M, 2, g]), tvals)
	
	p = plot(
		# title="Infectious (all groups)", 
		xlabel="",
		legend=:topright,
		frame=:box,
		yscale = :log10, 
		ylim = (I0/2, 2e-2),
	)
	plot!(tvals, X_wt, label="Wild type")
	plot!(tvals, X_m, label="Variant")
	
	p
end

# ╔═╡ 882379b0-76c3-4d74-b1f3-167f5c3c6900
pS_many = let
	g = :S
	tvals = range(0, T, length=400)
	
	X_wt = map(t -> mean(sol_many[t, 1:M, 1, g]), tvals)
	X_m = map(t -> mean(sol_many[t, 1:M, 2, g]), tvals)
	
	p = plot(
		# title="Infectious (all groups)", 
		xlabel="",
		legend=:topright,
		frame=:box,
		# yscale = :log10, 
		ylim = S_yrange,
	)
	plot!(tvals, X_wt, label="Wild type")
	plot!(tvals, X_m, label="Variant")
	hline!([δ/α], line = (10, .6, :gray), label="")
	annotate!(0.92*T, δ/α / 1.06, ("δ/α", 40))
	p
end

# ╔═╡ e6497af1-d7cf-4ad9-8b17-ef12b0853488
state_init_many_varC, sol_init_many_varC, sol_many_varC = let
	state = PSS.SIRState(; regions=many_groups, parameters=params_many_varC)
	# simulating for  time T without the variant to get initial state
	sol_init = PSS.simulate(state, (0, T))
	# introduce the variant
	state_init = PSS.set_infected(sol_init(T), 2, I0)
	# re-simulate for time T
	sol = PSS.simulate(state_init, (0, T))

	state_init, sol_init, sol
end;

# ╔═╡ 54244052-2b95-4e72-ac24-b5923d3081fa
p_freq_many_varC = let
	# freq. plot
	tvals = range(0, T, length=1000)
	f_R1 = PSS.frequency(sol_many_varC, tvals, 1, 2)
	f_R = PSS.frequency(sol_many_varC, tvals, 2)

	p = plot(
		legend = :bottomright,
		xlabel = "Time",
		ylabel = "",
		ylim = (-0.01, 1.01),
	)

	plot!(tvals, f_R, line=(:blue), label="")
		
	hline!([β]; label="", line=(:gray, 10, .5))

	annotate!(T*.95, β+.1, ("β", 40))
	annotate!(T*.15, .17, ("M=$M immune groups", 30, :left, "Helvetica Bold"))
	annotate!(T*.15, .07, ("Slow mixing", 30, :left, "Helvetica Bold"))
	
	p
end

# ╔═╡ feb7c4cd-d7b7-413b-8300-8c86283731ca
pI_many_varC = let
	g = :I
	tvals = range(0, T, length=400)
	
	X_wt = map(t -> mean(sol_many_varC[t, 1:M, 1, g]), tvals)
	X_m = map(t -> mean(sol_many_varC[t, 1:M, 2, g]), tvals)
	
	p = plot(
		# title="Infectious (all groups)", 
		xlabel="Time",
		legend=:topright,
		frame=:box,
		yscale = :log10, 
		ylim = (I0/2, 2e-2),
	)
	plot!(tvals, X_wt, label="Wild type")
	plot!(tvals, X_m, label="Variant")
	
	p
end

# ╔═╡ d7832f42-eb72-4483-9874-002af03f8773
pS_many_varC = let
	g = :S
	tvals = range(0, T, length=400)
	
	X_wt = map(t -> mean(sol_many_varC[t, 1:M, 1, g]), tvals)
	X_m = map(t -> mean(sol_many_varC[t, 1:M, 2, g]), tvals)
	
	p = plot(
		# title="Infectious (all groups)", 
		xlabel="Time",
		legend=:topright,
		frame=:box,
		# yscale = :log10, 
		ylim = S_yrange,
	)
	plot!(tvals, X_wt, label="Wild type")
	plot!(tvals, X_m, label="Variant")
	hline!([δ/α], line = (10, .6, :gray), label="")
	annotate!(0.92*T, δ/α / 1.06, ("δ/α", 40))
	p
end

# ╔═╡ f617e5fd-1ecf-460c-89c0-44f4405cfaca
state_init_one, sol_init_one, sol_one = let
	state = PSS.SIRState(; regions=one_group, parameters=params_one)
	# simulating for  time T without the variant to get initial state
	sol_init = PSS.simulate(state, (0, T))
	# introduce the variant
	state_init = PSS.set_infected(sol_init(T), 2, I0)
	# re-simulate for time T
	sol = PSS.simulate(state_init, (0, T))

	state_init, sol_init, sol
end;

# ╔═╡ cd5131c3-2a6f-46e9-9cdc-193573e516f8
p_freq_one = let
	# freq. plot
	tvals = range(0, T, length=100)
	f_R = PSS.frequency(sol_one, tvals, 2)

	p = plot(
		legend = :bottomright,
		# xlabel = "Time",
		ylabel = "",
		title = "Variant frequency",
		ylim = (-0.01, 1.01),
	)

	plot!(tvals, f_R, line=(:blue), label="")
	hline!([β]; label="", line=(:gray, 10, .5))

	annotate!(T*.95, β+.1, ("β", 40))
	annotate!(T*.15, .07, ("One immune group", 30, :left, "Helvetica Bold"))
	# plot logistic growth for comparison
	plot!(tvals, logistic_growth(
		tvals, growth_rate_one, f_R[1]);
		label = "", line=(:black, :dash, 2),
	)
	
	p
end

# ╔═╡ 1d6ff2eb-92dd-4588-ba4c-6114231c8ef4
pI_one = let
	g = :I
	tvals = range(0, T, length=400)
	
	X_wt = map(t -> sol_one[t, 1, 1, g], tvals)
	X_m = map(t -> sol_one[t, 1, 2, g], tvals)
	
	p = plot(
		title="Infectious", 
		# xlabel="Time",
		legend=:topright,
		frame=:box,
		yscale = :log10, 
		ylim = (I0/2, 2e-2),
	)
	plot!(tvals, X_wt, label="Wild type")
	plot!(tvals, X_m, label="Variant")
	
	p
end

# ╔═╡ 74daa07e-62f3-4937-af99-a726933bd647
pS_one = let
	g = :S
	tvals = range(0, T, length=400)
	
	X_wt = map(t -> sol_one[t, 1, 1, g], tvals)
	X_m = map(t -> sol_one[t, 1, 2, g], tvals)
	
	p = plot(
		title="Susceptibles", 
		# xlabel="Time",
		legend=:topright,
		frame=:box,
		# yscale = :log10, 
		ylim = S_yrange,
	)
	plot!(tvals, X_wt, label="Wild type")
	plot!(tvals, X_m, label="Variant")
	hline!([δ/α], line = (10, .6, :gray), label="")
	annotate!(0.92*T, δ/α / 1.06, ("δ/α", 40))
	p
end

# ╔═╡ dd51a9d9-7817-49d1-ae94-c69f4aac1c47
let
	panel = plot(
		pI_one, pS_one, p_freq_one, 
		pI_many, pS_many, p_freq_many,
		pI_many_varC, pS_many_varC, p_freq_many_varC,;
		layout = grid(3, 3),
		size = (2400, 2000),
		bottom_margin = 10mm, 
		right_margin = 10mm,
		dpi = 300, 
	)
	savefig("../figures/variant_invasion_SIR.png")
	savefig(joinpath(
		homedir(), 
		"Documents/BaleLabo/Notes/ExpiringFitness/figures/variant_invasion_SIR.png"
	))
	panel
end

# ╔═╡ 6d415fbc-b605-4a7c-b17f-db6480373e70
Cs

# ╔═╡ Cell order:
# ╠═4982e04e-6501-11ed-132e-ad3e395233da
# ╠═c2907ee7-3c4e-4b6b-946b-85971c5abbff
# ╟─6ab4e0ec-8ff1-4415-bcb8-5e5fbaf45bfc
# ╠═5e0b0485-e3dc-401a-818f-bcd8924995b5
# ╠═f905453d-2b68-40fd-9950-a7165933e73f
# ╟─db663438-8523-47ac-a783-88ff59a7ab1d
# ╟─c6c407d9-fea7-4e79-97ff-6acb88f4356d
# ╠═9c929838-c809-42c2-81b3-546a9de52e94
# ╠═98f8262c-28e4-4be1-91e9-6ac96bbb72a8
# ╠═f431ba40-9742-4898-92f0-33c3b8e12a1c
# ╠═e38f872f-f9ec-45cf-953c-42a18aa0bab3
# ╠═c2a159f4-cd25-4525-a270-2b3065a0dc0e
# ╠═943ee06f-d4e2-40fe-9764-6c0a0b7ef0b0
# ╟─a7fc518a-f6d3-4e93-bac5-02d946c768f8
# ╠═c4da46f5-8178-4ae5-b865-db025a365da0
# ╠═732be632-4b69-4434-bf51-33a219632fe2
# ╟─81cbea6e-9a97-4e7e-9528-f0830c7ae58f
# ╠═78b03821-37d4-4a9d-8ff0-349e12f66686
# ╟─b7dd719a-118e-4809-84da-609735534835
# ╠═a62f2142-7d76-4a15-9566-54115aa98e08
# ╟─51c77696-d4ae-466e-a19a-cf9278c7d3a5
# ╟─559e0978-3bb1-407f-b27b-4e31dc9df04e
# ╠═fd0d5424-074d-4fee-a41c-a80e61591348
# ╠═e6497af1-d7cf-4ad9-8b17-ef12b0853488
# ╠═f617e5fd-1ecf-460c-89c0-44f4405cfaca
# ╟─64c6f22b-9717-46c4-93cf-e1eafd107931
# ╟─1d3084e3-a18e-4850-b424-a7a32f4780b3
# ╟─7aa5f567-8af0-4b40-8669-f959da124826
# ╟─ffebafd7-5abf-461f-af8d-030448f77f44
# ╠═0e026c1f-a6c9-49d5-8cf2-23d0b33a217f
# ╠═29fa649f-7400-4099-8163-75a3d12362c9
# ╟─6c3a3b73-34b2-4643-9170-ed644e363725
# ╠═6d415fbc-b605-4a7c-b17f-db6480373e70
# ╟─dd51a9d9-7817-49d1-ae94-c69f4aac1c47
# ╟─d65daa41-6477-4eef-8193-364ee4cc21f6
# ╠═cd5131c3-2a6f-46e9-9cdc-193573e516f8
# ╠═d1910dd1-3e5e-4d1c-92cb-529b9a00c356
# ╠═d33785f7-96b6-4c26-92a7-82314b8114e1
# ╠═54244052-2b95-4e72-ac24-b5923d3081fa
# ╠═d2b8e952-37e1-4332-a321-706bd589a62d
# ╠═8f365998-cb7e-4000-a5c5-516a2b1ac814
# ╟─38e8868e-43ff-4f8b-88fb-07930638c0dc
# ╟─1d6ff2eb-92dd-4588-ba4c-6114231c8ef4
# ╟─2f8e4f1f-2da0-4bfa-b63c-069082180b6d
# ╟─feb7c4cd-d7b7-413b-8300-8c86283731ca
# ╟─d4d20d93-d54b-4fb8-a334-dde0a6ef5684
# ╠═cad2211c-3f9f-43e7-8098-17c07dc4c12f
# ╠═74daa07e-62f3-4937-af99-a726933bd647
# ╠═882379b0-76c3-4d74-b1f3-167f5c3c6900
# ╟─d7832f42-eb72-4483-9874-002af03f8773
# ╟─4640dbe6-89b5-4590-ac21-914226ca2097
# ╠═c79fe879-a32a-4391-8df1-5b9d9eeb75f5
