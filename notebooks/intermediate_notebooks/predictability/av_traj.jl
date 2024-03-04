### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 5c00618c-89c3-11ee-0991-95347f29c98d
begin
	using Pkg; Pkg.activate("../../../")
	using CSV
	using Chain
	using DataFrames
	using Distributions
	using Measures
	using Plots
	using StatsBase
end

# ╔═╡ 0f10333e-bd5c-4faa-97af-e2f56518502f
readdir("./trajectories_expiring_fitness/")

# ╔═╡ 83e60b11-e4e0-4cc7-aecb-43d8164765af


# ╔═╡ 993b804a-e640-4890-bb39-ed00c3b528ec


# ╔═╡ 91fcf1a6-0c6b-4965-b2b8-1cd50348384e


# ╔═╡ a2128c78-92aa-4836-b04c-a9f48de7f84e
S_to_vec(S) = @chain split(S, r",| |\[|\]") filter!(!isempty, _) parse.(Float64, _)

# ╔═╡ 5aa85583-e571-4d25-92d2-9b61a99fbf11
data_fixed = let
	df = DataFrame(CSV.File(
		"trajectories_fixed_fitness/data_avtraj.csv"
	))
	select!(
		df, Not(:mt, :mf), 
		[:mt, :mf] => ByRow((x1, x2) -> (S_to_vec(x1), S_to_vec(x2))) => [:mt, :mf]
	)
	df.α = zeros(Float64, size(df, 1))
	df
end;

# ╔═╡ ce2641d7-aa59-4ac0-b07d-bee9630b76ce
data_ef = let
	df = DataFrame(CSV.File(
		"trajectories_expiring_fitness/data_avtraj.csv"
	))
	select!(
		df, Not(:mt, :mf), 
		[:mt, :mf] => ByRow((x1, x2) -> (S_to_vec(x1), S_to_vec(x2))) => [:mt, :mf]
	)
end;

# ╔═╡ 2f0540f7-3192-4a9f-909c-5a99865f835e
# data = vcat(data_fixed, data_ef);
data = data_ef;

# ╔═╡ a571598b-9b62-4e0f-8347-f9aa88c52197
begin
	αvals = data.α |> sort |> unique
	ρvals = data.ρ |> sort |> unique
	s0 = data.s[1]
end

# ╔═╡ 785aaa23-cd8f-4658-9236-3d842b866bd3
let
	ρ = ρvals[1]
	k = 1
	
	dat = sort(data[data.ρ .== ρ, :], [:α])
	pal = palette(:bluesreds, length(αvals)+1)
	
	p = plot(
		xlabel = "time",
		ylabel = "frequency",
		frame=:box,
		title = "",
		legend = (k == 1 ? :topleft : false),
		xlim = (-100, 1000),
		ylim = (-0.01, 1.01),
	)
	
	annotate!(.8, .05, text("ρ/s = $(ρ/s0)", 24))
	
	for (i, r) in enumerate(eachrow(dat))
		plot!(
			p, r.mt, r.mf;
			label="α/s=$(round(r.α/r.s; sigdigits=2))",
			linewidth = 2,
			color = pal[i],
		)
	end

	p
end

# ╔═╡ Cell order:
# ╠═5c00618c-89c3-11ee-0991-95347f29c98d
# ╠═0f10333e-bd5c-4faa-97af-e2f56518502f
# ╠═5aa85583-e571-4d25-92d2-9b61a99fbf11
# ╠═ce2641d7-aa59-4ac0-b07d-bee9630b76ce
# ╠═83e60b11-e4e0-4cc7-aecb-43d8164765af
# ╠═2f0540f7-3192-4a9f-909c-5a99865f835e
# ╠═785aaa23-cd8f-4658-9236-3d842b866bd3
# ╠═993b804a-e640-4890-bb39-ed00c3b528ec
# ╠═91fcf1a6-0c6b-4965-b2b8-1cd50348384e
# ╠═a571598b-9b62-4e0f-8347-f9aa88c52197
# ╠═a2128c78-92aa-4836-b04c-a9f48de7f84e
