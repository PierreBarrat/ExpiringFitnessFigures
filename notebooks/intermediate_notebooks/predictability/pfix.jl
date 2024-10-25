### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ 7697c23c-887a-11ee-2ea0-618c103b9196
begin
	using Pkg; Pkg.activate("../../../")
	using CSV
	using Chain
	using DataFrames
	using Distributions
	using Measures
	using QuadGK
	using Plots
	using StatsBase
end

# ╔═╡ 0baab00c-7447-424c-a8d0-9753ad1035a8
S_to_vec(S) = @chain split(S, r",| |\[|\]") filter!(!isempty, _) parse.(Float64, _)

# ╔═╡ d534d5fc-a321-40d3-93e2-b6298bed3f9c
data_fixed = let
	df = DataFrame(CSV.File(
		"trajectories_fixed_fitness//data_pfix.csv"
	))
	select!(
		df, Not(:pfix, :f), 
		[:f, :pfix] => ByRow((x1, x2) -> (S_to_vec(x1), S_to_vec(x2))) => [:f, :pfix]
	)
	df.α = zeros(Float64, size(df, 1))
	df
end;

# ╔═╡ c41ec997-4eb5-4566-8b50-94c6ea7e7b82
data_ef = let
	df = DataFrame(CSV.File(
		"trajectories_expiring_fitness/data_pfix.csv"
	))
	select!(
		df, Not(:pfix, :f), 
		[:f, :pfix] => ByRow((x1, x2) -> (S_to_vec(x1), S_to_vec(x2))) => [:f, :pfix]
	)
end;

# ╔═╡ edd312e5-ede5-4217-a9f1-f3038072da2d
# data = vcat(data_ef, data_fixed);
data = data_ef;

# ╔═╡ 2f8dc8a6-0a02-4983-8b14-0db6af50d37a
data[data.α./data.s .== 3.0, :][1,:].pfix[2]

# ╔═╡ 847f6080-fd36-4baa-b978-7835420ee0ae
begin
	αvals = data.α |> sort |> unique
	ρvals = data.ρ |> sort |> unique
	s0 = data.s[1]
end

# ╔═╡ 0583c8d1-b309-4e11-b437-d7304338834f
plts = map(enumerate(ρvals)) do (k, ρ)
	dat = sort(data[data.ρ .== ρ, :], [:α])
	pal = palette(:bluesreds, length(αvals))
	
	p = plot(
		xlabel = "frequency",
		ylabel = "",
		frame=:box,
		title = "",
		legend = (k == 1 ? :bottomleft : false),
	)
	annotate!(.8, .05, text("ρ/s = $(ρ/s0)", 24))
	for (i, r) in enumerate(eachrow(dat))
		plot!(
			p, r.f, r.pfix;
			label="α/s=$(round(r.α/r.s; sigdigits=2))", marker=:cross,
			linewidth = 5,
			markersize = 10,
			color = pal[i],
		)
		# plot!(p, r.f, r.f .+ r.mβ*(1 .- r.f), label="", color=pal[i], linestyle=:dash)
	end
	plot!([0,1], [0,1], line=(:black, :dash), label="")
end

# ╔═╡ 3c06ee20-577e-4fbb-a061-478fc1250672
plot(
	plts..., 
	layout = grid(1, length(ρvals)), 
	size = (length(ρvals,)*600, 600),
	# plot_title=text("Probability of fixation", 20),
	bottommargin = 20mm,
	legendfontsize = 18,
	tickfontsize = 18,
	annotationfontsize=28, 
)

# ╔═╡ 93261d1d-3674-4fca-a6b1-ee857e6b3705
plts_ef = map(αvals) do α
	dat = sort(data_ef[data_ef.α .== α, :], [:ρ])

	p = plot(
		xlabel = "f",
		ylabel = "pfix",
	)

	for r in eachrow(dat)
		scatter!(p, r.f, r.pfix, label="s/ρ=$(s0/r.ρ)")
	end
	plot!([0,1], [0,1], line=(:black, :dash), label="")
end

# ╔═╡ e236a2ae-05c8-48d1-9ed0-d6a22e4e5884
plts_ef[3]

# ╔═╡ a9ddf8e5-11e5-4ba0-874b-4b6a59b56639
function estimate_pfix(x, s0, α)
	P(β) = (1-β)^(α/s0 - 1)
	eps = 1e-4
	Z = quadgk(β -> P(β), eps, 1-eps)[1]

	return x * quadgk(β -> P(β)/Z, eps, x)[1] + quadgk(β -> β*P(β)/Z, x, 1-eps)[1]
end

# ╔═╡ 0c9b99eb-2e73-432d-aef1-0a1db2b1ad64
estimate_pfix(.5, s0, .09)

# ╔═╡ Cell order:
# ╠═7697c23c-887a-11ee-2ea0-618c103b9196
# ╠═d534d5fc-a321-40d3-93e2-b6298bed3f9c
# ╠═c41ec997-4eb5-4566-8b50-94c6ea7e7b82
# ╠═edd312e5-ede5-4217-a9f1-f3038072da2d
# ╠═0583c8d1-b309-4e11-b437-d7304338834f
# ╠═2f8dc8a6-0a02-4983-8b14-0db6af50d37a
# ╠═3c06ee20-577e-4fbb-a061-478fc1250672
# ╠═93261d1d-3674-4fca-a6b1-ee857e6b3705
# ╠═e236a2ae-05c8-48d1-9ed0-d6a22e4e5884
# ╠═847f6080-fd36-4baa-b978-7835420ee0ae
# ╠═0baab00c-7447-424c-a8d0-9753ad1035a8
# ╠═a9ddf8e5-11e5-4ba0-874b-4b6a59b56639
# ╠═0c9b99eb-2e73-432d-aef1-0a1db2b1ad64
