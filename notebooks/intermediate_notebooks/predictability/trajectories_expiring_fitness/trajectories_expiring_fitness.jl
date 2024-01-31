### A Pluto.jl notebook ###
# v0.19.32

using Markdown
using InteractiveUtils

# ╔═╡ 10b7789c-e56f-4b71-8e38-312af11bb751
begin
	using Pkg; Pkg.activate("../../../../")
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

# ╔═╡ 904296da-d2ee-4d1c-8962-4144b013745d


# ╔═╡ d8666482-d8a4-43f4-ad21-81f6e972ab15
savedir = "data_" * notebook_name * "/"

# ╔═╡ 4a8fa33e-f996-4cf9-bd95-05de2bbd11f5
mkpath(savedir)

# ╔═╡ 63a63a7c-fa1c-4db8-8ce8-b7651d589dc0


# ╔═╡ 0ed067a4-3780-432b-bd54-9d52a22f2bbf
begin
    N = 50_000
    L = 400
    selected_sites = 1:L
    neutral_sites = 0:-1
    μ = 0
    sweep_times = :random
    md"**Constants**"
end

# ╔═╡ 056a5776-cf9f-4c90-8435-534ab7252f3c
s0 = 0.03

# ╔═╡ aa039c92-7a8b-4acc-8100-0f808d8fc6d1
begin
    ρvals = [1/20, 1/5, 1] * s0
    αvals = [0., 1/3, 1, 3] * s0
end

# ╔═╡ 8a82ed87-ca19-4afd-b6ed-090d1bdc20a7
md"## Writing to files"

# ╔═╡ 3ad99196-c1f4-4153-8c37-9b5f8db3bd39
md"# Functions"

# ╔═╡ e2097545-78dc-423a-81fc-bff679176875
get_Ne(ρ, β2) = (sweep_times == :periodic ? 1/(2*ρ)/β2 : 1/ρ/β2)

# ╔═╡ 1110d7fa-9c96-4b40-ae60-2396a08d4db0
parameters = map(Iterators.product(αvals, ρvals)) do (α, ρ)
    mβ = s0 / (s0 + α)
    β2 = s0/(s0 + α) * 2*s0 / (2*s0 + α) # for exponential distribution of s at fixed α
    Ne = get_Ne(ρ, β2)
    T = min(round(Int, 2_000*Ne), 1_000_000)
    Δt = 10
	(α=α, ρ=ρ, mβ = mβ, β2 = β2, T = T, Δt = Δt)
end |> x -> vcat(x...)

# ╔═╡ bbb24ac7-d514-4343-8c07-5f017ce45e0b
function simulate(
    α, ρ; Δt = 1, T = 10_000,
)
    # setting parameters
    fitness_distribution = Exponential(s0/2)
    switchgen = ρ > 0 ? round(Int, 1/ρ) : Inf

    # initial population
    H = zeros(Float64, L)
    ϕ = ExpiringFitness(; L=L, H, α)
    μ = 0.
    pop = Pop(ϕ; N, L, μ, sampling_method=:multinomial)

    # callbacks
    cb = (
        f1 = WF.frequencies,
        # polymorphism = pop -> WF.diversity(
        #     pop; method = :polymorphism, positions = neutral_sites,
        # ),
        varpos_strict = pop -> WF.diversity(
            pop; method = :variable_positions, variable=1/N, positions=selected_sites,
        ),
        varpos_5 = pop -> WF.diversity(
            pop; method = :variable_positions, variable=.05, positions=selected_sites,
        ),
        n_genotypes = pop -> length(pop.genotypes),
    )


    β2 = s0/(s0 + α) * 2*s0 / (2*s0 + α)
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
            varpos_5 = x.varpos_5,
            varpos_strict = x.varpos_strict,
            n_genotypes = x.n_genotypes,
        )
    end

    f1 = zeros(Float64, length(cb_vals), 2*L)
    for (i, x) in enumerate(cb_vals)
        f1[i, :] .= x.f1
    end
    trajectories = get_trajectories(f1; time_scaling=Δt) |> FrequencyTrajectories.trajectories_to_dataframe

    return trajectories, diversity
end

# ╔═╡ cb6caaf7-1088-49ff-afcd-b33cba0107ee
trajectories, diversity = let
    tjs = Dict()
    div = Dict()
    for (i, p) in enumerate(parameters)
        @info p i/length(parameters)
        @time Ts, d = simulate(p.α, p.ρ; Δt=p.Δt, T=p.T)
        tjs[p] = Ts

        df = DataFrame()
        for x in d
            push!(df, x)
        end
        div[p] = df
    end
    tjs, div
end;


# ╔═╡ 96f68927-a67c-46be-98e4-e0343db40c3d
filenames = map(enumerate(collect(keys(diversity)))) do (idx, p)
	idx => (
		α = p.α,
		ρ = p.ρ,
		s = s0,
		mβ = p.mβ,
		β2 = p.β2,
        Ne = get_Ne(p.ρ, p.β2),
        sweep_times = sweep_times,
		Δt = p.Δt,
        L=L,
		μ=μ,
		N=N,
        trajectory_file = "trajectory_$(idx).csv",
		diversity_file = "diversity_$(idx).csv",
		idx = idx,
		key = p,
	)
end

# ╔═╡ ef370192-358c-4610-9b8a-cc2c2cfefd16
open(savedir * "/files.json", "w") do io
	JSON3.pretty(io, Dict(filenames))
end

# ╔═╡ 90916994-59e1-459e-afd5-1d9fe8d4250a
for (i, v) in filenames
    open(savedir * v.diversity_file, "w") do io
        CSV.write(io, diversity[v.key])
    end
    open(savedir * v.trajectory_file, "w") do io
        CSV.write(io, trajectories[v.key])
    end
end

# ╔═╡ 344605b7-7735-49e2-b965-37dac463ab85
md"# Tests"

# ╔═╡ 97690e71-94ea-42a8-b0a8-b60a73d3ac3c
isfinite

# ╔═╡ Cell order:
# ╠═10b7789c-e56f-4b71-8e38-312af11bb751
# ╠═19759f0e-ced3-11ed-2b40-8b01992454e2
# ╠═904296da-d2ee-4d1c-8962-4144b013745d
# ╠═4a8fa33e-f996-4cf9-bd95-05de2bbd11f5
# ╠═d8666482-d8a4-43f4-ad21-81f6e972ab15
# ╠═63a63a7c-fa1c-4db8-8ce8-b7651d589dc0
# ╠═0ed067a4-3780-432b-bd54-9d52a22f2bbf
# ╠═056a5776-cf9f-4c90-8435-534ab7252f3c
# ╠═aa039c92-7a8b-4acc-8100-0f808d8fc6d1
# ╠═1110d7fa-9c96-4b40-ae60-2396a08d4db0
# ╠═bbb24ac7-d514-4343-8c07-5f017ce45e0b
# ╠═cb6caaf7-1088-49ff-afcd-b33cba0107ee
# ╠═8a82ed87-ca19-4afd-b6ed-090d1bdc20a7
# ╠═96f68927-a67c-46be-98e4-e0343db40c3d
# ╠═ef370192-358c-4610-9b8a-cc2c2cfefd16
# ╠═90916994-59e1-459e-afd5-1d9fe8d4250a
# ╟─3ad99196-c1f4-4153-8c37-9b5f8db3bd39
# ╠═e2097545-78dc-423a-81fc-bff679176875
# ╠═344605b7-7735-49e2-b965-37dac463ab85
# ╠═97690e71-94ea-42a8-b0a8-b60a73d3ac3c
