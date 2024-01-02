### A Pluto.jl notebook ###
# v0.19.36

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

# ╔═╡ 4c339e0c-98a8-11ed-06d4-6d227a4cfad5
begin
	using Revise
	using Pkg; Pkg.activate("../ArticleFigures/")
	using DiscreteRandomWalk1D
	using FrequencyTrajectories
	using Measures
	using PlutoUI
	using Plots
	using StatsBase
end

# ╔═╡ 1656f29a-b8d5-43a9-829e-799a43144a79
simulate = @bind clicked Button("Simulate")

# ╔═╡ e4780102-964f-4f19-bd1c-ff58811fdcce
eval_include = @bind click2 CheckBox()

# ╔═╡ 47b5c308-9ad6-49b5-871d-a71c000a2a77
begin
	x0 = .5
	β = .3
	ρ = .1
end

# ╔═╡ 6aebac0f-1ee0-43c8-881f-10bf8eb7c5d3
begin
	nreps = 3
	dx = [-0.05, 0, 0.05]
	dt = [-0.05, 0, 0.05]
end

# ╔═╡ 2a34c468-3c82-4290-a490-304a2169106c
T = 3*Int(round(1/β/β/ρ)) # simulation time

# ╔═╡ df7e5f19-9b40-478a-a86e-18257df25d81
timevals, X = let 
	clicked
	X = Matrix{Float64}(undef, nreps, T+1)
	timevals = nothing
	for r in 1:nreps
		traj = RW.trajectory!(RW.EFRW(β; x=x0 + dx[r], ρ), T)
		X[r, :] .= traj.x
		timevals = traj.t
	end
	timevals, X
end

# ╔═╡ 4b823bf4-c27c-4337-8e70-b8125c764e11
save = @bind clicked_save CheckBox()

# ╔═╡ ba170ee5-3fa4-45c5-a2bf-513c78a1c367
let
	sv = homedir() * "/Documents/BaleLabo/Notes/ExpiringFitness/figures/"
	clicked_save && savefig(sv * "example_random_walk.png")
end

# ╔═╡ 0a2f050b-fb7a-478c-8cda-e41f02cb03d3
simulate

# ╔═╡ e8bed50c-4a60-45cb-a229-7b4a9339c416
md"## More repeats"

# ╔═╡ 44ad1f71-667d-4bba-a48e-46bf7160c573
timevals_long, X_long = let 
	nreps = 1_00
	T = 3*Int(round(1/β/β/ρ))
	X = Matrix{Float64}(undef, nreps, T+1)
	timevals = nothing
	for r in 1:nreps
		traj = RW.trajectory!(RW.EFRW(β; x=x0, ρ), T)
		X[r, :] .= traj.x
		timevals = traj.t
	end
	timevals, X
end

# ╔═╡ 6f5338c2-d64f-4e4e-97ce-8d3bee005814
mean(X_long, dims=1)

# ╔═╡ 85f077a3-70cd-42d8-97da-2dd5307c8201
timevals_long

# ╔═╡ c2051d10-5d95-4251-a62b-dc81b6fadd6d
md"## Helper functions"

# ╔═╡ ab29f9d5-4f99-4583-908a-d11b612bc9bd
function final_state_color(x::AbstractVector)
	return if x[end] > 0.95
		:red
	elseif x[end] < 0.05
		:blue
	else
		:gray
	end
end

# ╔═╡ 11958921-7425-4ec8-9a95-124171d69771
let
	p = plot()
	for freq in eachrow(X_long)
		color = final_state_color(freq)
		plot!(timevals_long, freq; color, label="", alpha = .1, linewidth=1)
	end

	plot!(timevals_long, vec(mean(X_long, dims=1)), line=(4, :black), label="")
	p
end

# ╔═╡ 55666413-027a-4b00-88bf-7bcc8299310a
function flat_interp(X, times)
	real_times = 0:length(X)
	X_interp = similar(times)
	
	pos = 1
	for (i,t) in enumerate(times)
		if t > real_times[pos+1]
			pos += 1
		end
		X_interp[i] = X[pos]
	end
	X_interp
end

# ╔═╡ 94226348-3a21-4aa8-9b7b-61c4922c17be
plt = let
	p = plot(
		xlabel = "Time",
		# xlim = (-0.2, 8),
	)
	tvals = range(extrema(timevals)..., length=1000)
	for (r, freq) in enumerate(eachrow(X))
		x = flat_interp(freq, tvals)
		plot!(tvals .+ dt[r], x, label="", linewidth=5)
	end
	p
end

# ╔═╡ Cell order:
# ╠═4c339e0c-98a8-11ed-06d4-6d227a4cfad5
# ╠═1656f29a-b8d5-43a9-829e-799a43144a79
# ╠═e4780102-964f-4f19-bd1c-ff58811fdcce
# ╠═47b5c308-9ad6-49b5-871d-a71c000a2a77
# ╠═6aebac0f-1ee0-43c8-881f-10bf8eb7c5d3
# ╠═2a34c468-3c82-4290-a490-304a2169106c
# ╠═df7e5f19-9b40-478a-a86e-18257df25d81
# ╠═94226348-3a21-4aa8-9b7b-61c4922c17be
# ╠═4b823bf4-c27c-4337-8e70-b8125c764e11
# ╠═ba170ee5-3fa4-45c5-a2bf-513c78a1c367
# ╠═0a2f050b-fb7a-478c-8cda-e41f02cb03d3
# ╠═e8bed50c-4a60-45cb-a229-7b4a9339c416
# ╠═44ad1f71-667d-4bba-a48e-46bf7160c573
# ╠═11958921-7425-4ec8-9a95-124171d69771
# ╠═6f5338c2-d64f-4e4e-97ce-8d3bee005814
# ╠═85f077a3-70cd-42d8-97da-2dd5307c8201
# ╠═c2051d10-5d95-4251-a62b-dc81b6fadd6d
# ╠═ab29f9d5-4f99-4583-908a-d11b612bc9bd
# ╠═55666413-027a-4b00-88bf-7bcc8299310a
