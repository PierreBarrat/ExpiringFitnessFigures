### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ 10b7789c-e56f-4b71-8e38-312af11bb751
begin
	using Revise
	using Pkg; Pkg.activate("../../../../")
	using Chain
	using CSV
	using DataFrames
	using Distributions
	using FrequencyTrajectories
	using JSON3
	using WrightFisher
	using StatsBase
end

# ╔═╡ 19759f0e-ced3-11ed-2b40-8b01992454e2
notebook_name = replace(@__FILE__, r"\.jl#==#.*" => "") |> basename

# ╔═╡ d8666482-d8a4-43f4-ad21-81f6e972ab15
savedir = "data_" * notebook_name * "/"

# ╔═╡ 4a8fa33e-f996-4cf9-bd95-05de2bbd11f5
mkpath(savedir)

# ╔═╡ 6cf21eab-ddb6-4c91-9e6a-d925d61702c2


# ╔═╡ 0ed067a4-3780-432b-bd54-9d52a22f2bbf
begin
	N = 1_000_000
	L = 150
	μ = 0
    f0 = 0.01
    max_freq = f0
	md"**Constants**"
end

# ╔═╡ afe7bf04-d998-46fc-afed-8e2cd21361d2
begin
    s0 = 0.03
	svals = [s0]
	ρvals = ρvals = [1/20, 1/5, 1] * s0
    Δtvals = [10]
    cfs_vals = [:random]
	md"**Variable parameters**"
end

# ╔═╡ 1110d7fa-9c96-4b40-ae60-2396a08d4db0
parameters = map(Iterators.product(cfs_vals, Δtvals, ρvals, svals)) do (cfs, Δt, ρ, s)
	(s=s, ρ=ρ, Δt=Δt, cfs=cfs)
end |> (x -> vcat(x...));

# ╔═╡ bbb24ac7-d514-4343-8c07-5f017ce45e0b
function simulate(s, ρ; Δt = 1, change_field_time)
	# setting parameters
	fitness_distribution = Exponential(s/2)

	T = max(1000/ρ, 1000/s)
	
	switchgen = round(Int, 1/ρ)

	# initial population
	H = rand(fitness_distribution, L)
	ϕ = AdditiveFitness(; L, H)
	pop = Pop(ϕ; N, L, μ)
	cb = (
		f1 = WF.frequencies,
		varpos_strict = pop -> WF.diversity(
			pop; method = :variable_positions, variable=1/N
		)/L,
		varpos_5 = pop -> WF.diversity(
			pop; method = :variable_positions, variable=0.05,
		)/L,
		n_genotypes = pop -> length(pop.genotypes),
	)
	
	cb_vals, switch_times = WF.Tools.evolve_sample!(
		pop, T, Δt, cb;
		fitness_distribution,
		switchgen, 
		change_init_field = true,
		change_field_time,
		max_freq,
        f0,
	)

	diversity = map(cb_vals) do x
		(varpos_strict = x.varpos_strict, varpos_5 = x.varpos_5)
	end

	f1 = zeros(Float64, length(cb_vals), 2*L)
	for (i, x) in enumerate(cb_vals)
		f1[i, :] .= x.f1
	end
	trajectories = get_trajectories(f1; time_scaling=Δt) |> FrequencyTrajectories.trajectories_to_dataframe

	return trajectories, diversity, switch_times
end

# ╔═╡ e61957e2-8c90-4e75-a52d-3b195ac31078
trajectories, diversity = let
	tjs = Dict()
	div = Dict()
	for (i, p) in enumerate(parameters)
		@info p i/length(parameters)
		@time T, d = simulate(p.s, p.ρ; Δt=p.Δt, change_field_time = p.cfs)
		tjs[p] = T
		
		df = DataFrame()
		for x in d
			push!(df, x)
		end
		div[p] = df
	end
	tjs, div
end;

# ╔═╡ 292f7538-7212-4b02-a0e8-28f1b6137602
filenames = map(enumerate(collect(keys(trajectories)))) do (idx, p)
	idx => (
		ρ = p.ρ,
		s = p.s,
		Δt = p.Δt,
		L = L,
		N = N,
        cfs = p.cfs,
		trajectory_file = "trajectory_$(idx).csv",
		diversity_file = "diversity_$(idx).csv",
		idx = idx,
		key = p,
	)
end

# ╔═╡ 1e8c8a5b-a767-46f2-98b9-62c6ce8dbb2d
open(savedir * "/files.json", "w") do io
	JSON3.pretty(io, Dict(filenames))
end

# ╔═╡ 8776e961-93e4-40fb-861b-b98abee83871
for (i, v) in filenames
	open(savedir * v.diversity_file, "w") do io
		CSV.write(io, diversity[v.key])
	end
	open(savedir * v.trajectory_file, "w") do io
		CSV.write(io, trajectories[v.key])
	end
end

# ╔═╡ Cell order:
# ╠═10b7789c-e56f-4b71-8e38-312af11bb751
# ╠═19759f0e-ced3-11ed-2b40-8b01992454e2
# ╠═4a8fa33e-f996-4cf9-bd95-05de2bbd11f5
# ╠═d8666482-d8a4-43f4-ad21-81f6e972ab15
# ╟─6cf21eab-ddb6-4c91-9e6a-d925d61702c2
# ╠═0ed067a4-3780-432b-bd54-9d52a22f2bbf
# ╠═afe7bf04-d998-46fc-afed-8e2cd21361d2
# ╠═1110d7fa-9c96-4b40-ae60-2396a08d4db0
# ╠═bbb24ac7-d514-4343-8c07-5f017ce45e0b
# ╠═e61957e2-8c90-4e75-a52d-3b195ac31078
# ╠═292f7538-7212-4b02-a0e8-28f1b6137602
# ╠═1e8c8a5b-a767-46f2-98b9-62c6ce8dbb2d
# ╠═8776e961-93e4-40fb-861b-b98abee83871
