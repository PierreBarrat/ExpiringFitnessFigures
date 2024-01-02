### A Pluto.jl notebook ###
# v0.19.30

using Markdown
using InteractiveUtils

# ╔═╡ 10b7789c-e56f-4b71-8e38-312af11bb751
begin
	using Pkg; Pkg.activate("../../../")
	using Chain
	using CSV
	using DataFrames
	using DelimitedFiles
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

# ╔═╡ 63a63a7c-fa1c-4db8-8ce8-b7651d589dc0


# ╔═╡ 056a5776-cf9f-4c90-8435-534ab7252f3c
begin
	mβvals = 0.3:0.3:0.9 |> collect
	β2vals = map(mβvals) do m
        [m^2 + (m - m^2)/3]
	end

	βvals = reverse([(mβ=m, β2=v) for (j, m) in enumerate(mβvals) for v in β2vals[j]])
	md"**Distributions of β**"
end

# ╔═╡ aa039c92-7a8b-4acc-8100-0f808d8fc6d1
ρvals = @chain range(log(1e-2), log(1e-1), length=5) exp.(_) collect reverse vcat(_, 3e-3)

# ╔═╡ 8a82ed87-ca19-4afd-b6ed-090d1bdc20a7
md"## Writing to files"

# ╔═╡ 3ad99196-c1f4-4153-8c37-9b5f8db3bd39
md"# Functions"

# ╔═╡ b7d9a8d0-476d-4d87-acd8-a06c7766198b
begin
	struct FitnessDistribution <: Sampleable{Univariate, Continuous}
		β :: Distribution
		α :: Float64
	end

	function FitnessDistribution(mβ, β2, α)
		@assert β2 < mβ < sqrt(β2) "Incorrect value <β>=$(mβ). The range of possible values is $(β2) < mβ < $(sqrt(β2))"

		return FitnessDistribution(get_β_distribution(mβ, β2), α)
	end
end

# ╔═╡ 880c9804-d0e3-4f71-a34c-fcef988e3953
begin
	import Base.rand
	function rand(rng::Distributions.AbstractRNG, x::FitnessDistribution)
		β = @chain rand(x.β) min(_, .999) max(_, .001) # for the log to be ok
		return -x.α * log(1-β)/2
	end
	flag = true
end

# ╔═╡ 0ed067a4-3780-432b-bd54-9d52a22f2bbf
begin
    flag ? 1+1 : 2
    N = 50_000
    L_sel = 250
    L_neutral = 100
    L = L_sel + L_neutral
    selected_sites = 1:L_sel
    neutral_sites = (L_sel+1):L

    μ = 0
    # μ_neutral = 0.1/N
    α = .3
    # Δt = 5
    sweep_times = :random

    get_Ne(ρ, β2) = (sweep_times == :periodic ? 1/(2*ρ)/β2 : 1/ρ/β2)
	get_μ_neutral(ρ, β2) = 0.01 / get_Ne(ρ, β2)
    md"**Constants**"
end

# ╔═╡ 1110d7fa-9c96-4b40-ae60-2396a08d4db0
parameters = map(Iterators.product(βvals, ρvals)) do (βv, ρ)
    Ne = get_Ne(ρ, βv.β2)
    T = max(round(Int, 100*Ne), 20_000)
    Δt = round(Int, T / 2e4)
	(mβ = βv.mβ, β2 = βv.β2, ρ=ρ, μ_neutral = get_μ_neutral(ρ, βv.β2), T = T, Δt = Δt)
end |> x -> vcat(x...)

# ╔═╡ cb6caaf7-1088-49ff-afcd-b33cba0107ee
diversity = let
	div = Dict()
	for (i, p) in enumerate(parameters)
		@info p i/length(parameters)
		@time diversity = simulate(p.mβ, p.β2, p.ρ, p.μ_neutral; p.Δt, p.T)
		# diversity
		df = DataFrame()
		for x in diversity
			push!(df, x)
		end
		div[p] = df
	end
	div
end;

# ╔═╡ bbb24ac7-d514-4343-8c07-5f017ce45e0b
function simulate(
    mβ, β2, ρ, μ_neutral;
    Δt = 1, T = 10_000,
)
    # setting parameters
    fitness_distribution = FitnessDistribution(mβ, β2, α)
    switchgen = ρ > 0 ? round(Int, 1/ρ) : Inf

    # initial population
    H = zeros(Float64, L)
    ϕ = ExpiringFitness(; L=L, H, α)
    μ_vector = vcat(μ*ones(L_sel), μ_neutral*ones(L_neutral))
    pop = Pop(ϕ; N, L, μ = μ_vector, sampling_method=:multinomial)

    # callbacks
    cb = (
        # f1 = WF.frequencies,
        polymorphism = pop -> WF.diversity(
            pop; method = :polymorphism, positions = neutral_sites,
        ),
        varpos_strict = pop -> WF.diversity(
            pop; method = :variable_positions, variable=1/N, positions=selected_sites,
        ),
        # varpos_5 = pop -> WF.diversity(
        #     pop; method = :variable_positions, variable=.05, positions=selected_sites,
        # ),
        n_genotypes = pop -> length(pop.genotypes),
    )

    @info "Simulating for $T generations (Ne = $(get_Ne(ρ, β2)))"
    cb_vals, switch_times = WF.Tools.evolve_sample!(
        pop, T, Δt, cb;
        fitness_distribution,
        switchgen,
        change_init_field=false,
        change_field_time = sweep_times,
        max_freq=0.,
        selected_positions=selected_sites,
        burnin = get_Ne(ρ, β2),
    )

    diversity = map(cb_vals) do x
        (
            t = x.t,
            varpos_strict = x.varpos_strict,
            polymorphism = x.polymorphism,
            n_genotypes = x.n_genotypes,
        )
    end

    return diversity
end

# ╔═╡ 96f68927-a67c-46be-98e4-e0343db40c3d
filenames = map(enumerate(collect(keys(diversity)))) do (idx, p)
	idx => (
		mβ = p.mβ,
		β2 = p.β2,
		ρ = p.ρ,
        Ne = get_Ne(p.ρ, p.β2),
        sweep_times = sweep_times,
		α = α,
		Δt = p.Δt,
		L_sel=L_sel,
		L_neutral=L_neutral,
        L=L,
		μ=μ,
		μ_neutral=p.μ_neutral,
		N=N,
		diversity_file = "diversity_$(idx).csv",
		idx = idx,
		key = p,
	)
end

# ╔═╡ 233eccf0-0a03-4350-8ea8-7da66b7fbdb3
map(filenames) do Y
	X = Y[2]
	(2*X.ρ * X.β2), X.μ_neutral
end

# ╔═╡ 9d5b8152-3abf-43d1-81c4-101b43417a87
filenames

# ╔═╡ ef370192-358c-4610-9b8a-cc2c2cfefd16
open(savedir * "/files.json", "w") do io
	JSON3.pretty(io, Dict(filenames))
end

# ╔═╡ 90916994-59e1-459e-afd5-1d9fe8d4250a
for (i, v) in filenames
	open(savedir * v.diversity_file, "w") do io
		CSV.write(io, diversity[v.key])
	end
	# open(savedir * v.trajectory_file, "w") do io
	# 	CSV.write(io, trajectories[v.key])
	# end
	# open(savedir * v.allele_freq_file, "w") do io
	# 	writedlm(io, allele_frequencies[v.key])
	# end
	# open(savedir * v.switch_info_file, "w") do io
	# 	CSV.write(io, switch_info[v.key])
	# end
end

# ╔═╡ 57ac0330-e128-4d9c-9b00-8144de833009
let
	B = FitnessDistribution(0.3, 0.12, 0.5)
	rand(B)
end

# ╔═╡ f77e8aa9-7843-4917-a745-0fa3e6611939
"""
	get_β_distribution(mβ, β2)

Return a Beta distribution with mean `mβ` and second moment `β2`.
"""
function get_β_distribution(mβ, β2)
	b = (mβ - β2)/(β2 - mβ^2)*(1-mβ)
	a = mβ/(1-mβ)*b
	return Beta(a,b)
end

# ╔═╡ caa421a1-a3d2-4e13-af68-27bc08229c5a
sfs(f1) = sum(x -> x*(1-x), f1[1:2:end])

# ╔═╡ 344605b7-7735-49e2-b965-37dac463ab85
md"# Tests"

# ╔═╡ 97690e71-94ea-42a8-b0a8-b60a73d3ac3c
isfinite

# ╔═╡ Cell order:
# ╠═10b7789c-e56f-4b71-8e38-312af11bb751
# ╠═19759f0e-ced3-11ed-2b40-8b01992454e2
# ╠═4a8fa33e-f996-4cf9-bd95-05de2bbd11f5
# ╠═d8666482-d8a4-43f4-ad21-81f6e972ab15
# ╠═63a63a7c-fa1c-4db8-8ce8-b7651d589dc0
# ╠═0ed067a4-3780-432b-bd54-9d52a22f2bbf
# ╠═056a5776-cf9f-4c90-8435-534ab7252f3c
# ╠═aa039c92-7a8b-4acc-8100-0f808d8fc6d1
# ╠═1110d7fa-9c96-4b40-ae60-2396a08d4db0
# ╠═bbb24ac7-d514-4343-8c07-5f017ce45e0b
# ╠═cb6caaf7-1088-49ff-afcd-b33cba0107ee
# ╠═233eccf0-0a03-4350-8ea8-7da66b7fbdb3
# ╠═9d5b8152-3abf-43d1-81c4-101b43417a87
# ╠═8a82ed87-ca19-4afd-b6ed-090d1bdc20a7
# ╠═96f68927-a67c-46be-98e4-e0343db40c3d
# ╠═ef370192-358c-4610-9b8a-cc2c2cfefd16
# ╠═90916994-59e1-459e-afd5-1d9fe8d4250a
# ╟─3ad99196-c1f4-4153-8c37-9b5f8db3bd39
# ╠═880c9804-d0e3-4f71-a34c-fcef988e3953
# ╠═57ac0330-e128-4d9c-9b00-8144de833009
# ╠═b7d9a8d0-476d-4d87-acd8-a06c7766198b
# ╠═f77e8aa9-7843-4917-a745-0fa3e6611939
# ╠═caa421a1-a3d2-4e13-af68-27bc08229c5a
# ╠═344605b7-7735-49e2-b965-37dac463ab85
# ╠═97690e71-94ea-42a8-b0a8-b60a73d3ac3c
