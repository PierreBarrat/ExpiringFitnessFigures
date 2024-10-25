### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ 4c339e0c-98a8-11ed-06d4-6d227a4cfad5
# begin
# 	using Revise
# 	using Pkg; Pkg.activate("../ArticleFigures/")
# 	using DiscreteRandomWalk1D
# 	using FrequencyTrajectories
# 	using Measures
# 	using PlutoUI
# 	using Plots
# 	using StatsBase
# end

function random_walk_example_plot(seed_val; nreps = 3)

    Random.seed!(seed_val)

    # ╔═╡ 47b5c308-9ad6-49b5-871d-a71c000a2a77
    begin
    	x0 = .5
    	β = .3
    	ρ = .1
    end

    # ╔═╡ 6aebac0f-1ee0-43c8-881f-10bf8eb7c5d3
    begin
    	# nreps = 3
        (dx, dt) = if nreps == 3
    	    ([-0.05, 0, 0.05], [-0.05, 0, 0.05])
        else
            d = range(0 - div(nreps, 2)*0.025, 0 + div(nreps,2) * 0.025, length=nreps)
            (d, d)
        end
    end

    # ╔═╡ 2a34c468-3c82-4290-a490-304a2169106c
    T = 3*Int(round(1/β/β/ρ)) # simulation time

    # ╔═╡ df7e5f19-9b40-478a-a86e-18257df25d81
    timevals, X = let
    	X = Matrix{Float64}(undef, nreps, T+1)
    	timevals = nothing
    	for r in 1:nreps
    		traj = RW.trajectory!(RW.EFRW(β; x=x0 + dx[r], ρ), T)
    		X[r, :] .= traj.x
    		timevals = traj.t
    	end
    	timevals, X
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

    return plt

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

# ╔═╡ Cell order:
# ╠═4c339e0c-98a8-11ed-06d4-6d227a4cfad5
# ╠═47b5c308-9ad6-49b5-871d-a71c000a2a77
# ╠═6aebac0f-1ee0-43c8-881f-10bf8eb7c5d3
# ╠═2a34c468-3c82-4290-a490-304a2169106c
# ╠═df7e5f19-9b40-478a-a86e-18257df25d81
# ╠═94226348-3a21-4aa8-9b7b-61c4922c17be
# ╠═55666413-027a-4b00-88bf-7bcc8299310a
