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


function pop_stack_plot(seed_val; μ=.25, λ=.25, dist = :uniform, N=7)


    Random.seed!(seed_val)

    # ╔═╡ 5094c0ff-bf9e-4476-b0f1-9b9941c1c452
    md"# Setup"

    # ╔═╡ 446f0274-6944-43a9-bad8-265aec6ada60
    sim_button = @bind clicked Button("Simulate")

    # ╔═╡ e1e91203-6a97-43f0-b031-036853f15ee1
    seed_picker = @bind val NumberField(1:100_000)

    # ╔═╡ 629df600-203e-4456-b1f0-1479af1a4199
    params = let
    	# Random.seed!(val)
    	clicked
    	M = 10
    	α = 3
    	γ = .01
    	c = 1/M
    	PSS.Parameters(; N, M, α, γ, c)
    end;

    # ╔═╡ 3279b428-b03e-48e5-97e8-06f8f8e4d80f
    I0 = 1e-5 # initial value of each new mutant

    # ╔═╡ 528bfbbe-0102-447d-a9fd-5e1f17e55238
    @unpack α, γ, δ, N, M = params;

    # ╔═╡ 30f888a7-d3f4-4089-9c1e-c39ada15444a
    begin
    	local teq = 4/γ
    	# Random.seed!(val)
    	itimes = @chain begin
    		# rand(Exponential(Δt), N-1)
    		ones(N)*teq
    		cumsum
    		_ .- teq # first introduction at 0
    		pushfirst!(_, -5*teq) # wt starts equilibrated
    	end
    	itimes[end] += teq # longer simulation for the last strain
    	"Times at which each new strain is introduced: $(itimes)"
    end

    # ╔═╡ e34bdc45-3bb1-4628-8e55-c40d54d4159a
    K = let
    	# Random.seed!(val)
        # μ, λ = .25, .25
        backward, forward = if dist == :expo
    	    ϵf_dist = Exponential(λ)
    	    ϵb_dist = Exponential(μ)
    	    backward = max.(.1, 1 .- rand(ϵb_dist, params.N))
    	    forward = max.(.1, 1 .- rand(ϵf_dist, params.N))
            backward, forward
        elseif dist == :uniform
            f_dist = Uniform(.3, .7)
            b_dist = Uniform(.3, .7)
            backward = rand(b_dist, params.N)
            forward = rand(f_dist, params.N)
            backward, forward
        end
        # first strain at 0.5
        # backward[2] = .5
        # forward[2] = .5
    	PSS.cross_immunity(backward, forward)
    end;

    # ╔═╡ b761969b-8fc1-4772-aeaf-f61c6525fcdc
    function make_region(K)
    	S = δ/α
    	I = 0
    	C = 0
    	R = 1 - S - I - C
    	viruses = []
    	for v in 1:N
    		push!(viruses, PSS.Virus(; S, I, C, R))
    	end
    	PSS.Region(; viruses, K)
    end

    # ╔═╡ f75d7c6b-ca6f-4dfe-aa00-d460e7776918
    special_region = make_region(K);

    # ╔═╡ 1237699b-1a27-449d-94c8-983013ce1907
    normal_regions = [make_region(ones(N,N)) for m in 2:M];

    # ╔═╡ 0defdd45-4fca-49e4-9216-fbb06d9e0968
    state_init = PSS.SIRState(;
    	regions=vcat(special_region, normal_regions...),
    	parameters=params
    );

    # ╔═╡ 406c2e67-20e4-4062-aa04-58113b7a06fa
    md"# Simulation"

    # ╔═╡ 261a3400-1415-4a67-99f7-4d6f2bc2d19d
    md"# Plots"

    # ╔═╡ fc603500-6420-4018-8b5f-a670ef128c5a
    sim_button

    # ╔═╡ 6571f401-fc4e-4d04-8f8a-83f69b579701
    md"# Utils"

    # ╔═╡ 36ca0fdb-d0d6-45ee-93f1-090f761e3191
    md"#### Extracting vectors of I, t, freq from array of solutions"

    # ╔═╡ 41cc2b47-cab0-4f4b-a243-f8adc58deb2d
    begin
    	function t_values(solutions)
    		tvals = Float64[]
    		for sol in solutions
    			append!(tvals, range(sol.tspan..., length=100))
    		end
    		return tvals
    	end
    	function compartment_values(solutions, n; g=:I)
    		X = Float64[]
    		tvals = Float64[]
    		for sol in solutions
    			ts = range(sol.tspan..., length=100)
    			append!(tvals, ts)
    			append!(X, sum(sol[ts, i, n, g] for i in 1:params.M))
    		end
    		return tvals, X
    	end
    	function compartment_values(solutions; g=:I)
    		Xs = [compartment_values(solutions, n)[2] for n in 1:length(solutions)]
    		tvals = t_values(solutions)
    		return tvals, Xs
    	end
    end

    # ╔═╡ dce2093c-1e9d-47db-8cd1-d47e7a5216cf
    """
    	simulate(state_init, tspan)

    Simulate starting from from `state_init` during `tspan`.
    """
    function simulate(state_init, tspan)
    	sol = PSS.simulate(state_init, tspan)
    	state_final = sol(tspan[end])
    	return sol, state_final
    end

    # ╔═╡ 9b70216c-13ec-4637-b536-ceb85c7346ab
    solutions, states = let
    	solutions = []
    	states = [state_init]
    	current_state = state_init

    	for n in 1:N
    		tspan = (itimes[n], itimes[n+1])
    		current_state = PSS.set_infected(current_state, n, I0)
    		sol, current_state = simulate(current_state, tspan)
    		push!(solutions, sol)
    		push!(states, current_state)
    	end
    	solutions, states
    end;

    # ╔═╡ 13b2937c-116e-4476-afe9-5a4fad087c25
    let
    	tvals, Is = compartment_values(solutions)

    	p = plot(
    		yscale=:log10,
    		ylim = (1e-7, 1),
    		title = "Infected",
    		xlabel = "Time",
    		legend = :bottomright
    	)

    	for n in 1:N
    		plot!(p, tvals, Is[n]; label = "Strain $n")
    	end
    	p
    end;

    # ╔═╡ 5e4835fb-3447-4ca4-831b-cabbc2db3528
    """
    	pick_next_strain(idx_wt, solutions)

    Decide whether the next strain comes from the wild-type (`true`) or the mutant (`false`) background. Based on results of simulation and on status of previous strains `idx_wt`.
    """
    function pick_next_strain(idx_wt, solutions)
    	V = length(idx_wt)
    	f_wt = @chain begin
    		solutions
    		getindex(V)
    		[PSS.frequency(_, [_.tspan[end]], v)[1] for v in 1:V]
    		getindex(idx_wt)
    		sum
    	end

    	return rand() < f_wt ? true : false

    	return f_wt
    end

    # ╔═╡ 136175fd-5bdf-46ae-8bf4-55d6b9b47f08
    # Choose which strain is mutant background, which is w.t. background
    viruses, idx_wt =  let
    	# Random.seed!(val)
    	idx_wt = [true, false]
    	for v in 3:N
    		push!(idx_wt, pick_next_strain(idx_wt, solutions))
    	end
    	# idx_wt = append!([true, false], rand(Bool, N-2))
    	idx_mut = .!idx_wt

    	# Some bookkeeping based on this
    	viruses = []
    	m = 1
    	w = 1
    	for n in 1:N
    		if idx_wt[n]
    			push!(viruses, (:wt, w))
    			w += 1
    		else
    			push!(viruses, (:mut, m))
    			m += 1
    		end
    	end

    	viruses, idx_wt
    end;

    # ╔═╡ 80fa524d-f0b1-470c-a7e8-f4857e33ee4a
    begin
    	idx_mut = .!idx_wt
    	N_wt = count(idx_wt)
    	N_mut = count(.!idx_wt)
    end

    # ╔═╡ 424d9f47-c380-4bf4-9dab-be01a4de4830
    # begin
    # 	pal_wt = palette(:blues, N_wt+1)
    # 	pal_mut = palette(:reds, N_mut+1)
    # end;

    # # ╔═╡ 941cb880-c300-4785-bf76-f816aea7bced
    # p_pop_stack = let
    # 	p = plot(
    # 		# yscale=:log10,
    # 		ylim = (1e-5, 1),
    # 		title = "Frequency of mutation A",
    # 		xlabel = "Time",
    # 		xlim = (-20, maximum(itimes)),
    # 	)


    # 	tvals, Is = compartment_values(solutions)
    # 	Z = sum(Is)
    # 	# Wt
    # 	Is_wt = [sum(Is[findall(idx_wt)[1:n]])./Z for n in 1:N_wt]
    # 	last = ones(length(tvals))
    # 	for (v, I) in enumerate(Is_wt)
    # 		plot!(
    # 			tvals, 1 .- I;
    # 			color = pal_wt[v], fillrange=last, linewidth=0, label="",
    # 		)
    # 		last = 1 .- I
    # 	end

    # 	# Mut
    # 	Is_mut = [sum(Is[findall(idx_mut)[1:n]])./Z for n in 1:N_mut]
    # 	last = zeros(length(tvals))
    # 	for (v, I) in enumerate(Is_mut)
    # 		next = I
    # 		plot!(
    # 			tvals, last;
    # 			color = pal_mut[v+1], linewidth=0, label="", fillrange = next,
    # 		)
    # 		last = next
    # 	end

    # 	plot!(tvals, Is_mut[end], line=(:black), label="")
    # 	p
    # end;

    # # ╔═╡ 34c7c53a-dec9-49ae-b2c7-e4a3bcc020c1
    # p_freq = let
    # 	tvals, Is = compartment_values(solutions)

    # 	p = plot(
    # 		yscale=:linear,
    # 		ylim = (1e-2, 1),
    # 		title = "Frequency of variants",
    # 		xlabel = "Time",
    # 		xlim = (-20, maximum(itimes)),
    # 		legendfontsize = 22
    # 	)

    # 	first_m = true
    # 	first_wt = true
    # 	for (n,v) in enumerate(viruses)
    # 		label = if first_m && v[1] != :wt
    # 			first_m = false
    # 			"Variant background"
    # 		elseif first_wt && v[1] == :wt
    # 			first_wt = false
    # 			"Wild-type background"
    # 		else
    # 			""
    # 		end
    # 		plot!(
    # 			p, tvals, Is[n] ./ sum(Is);
    # 			color = v[1]==:wt ? pal_wt[v[2]] : pal_mut[v[2]],
    # 			label, linewidth = 5,
    # 		)
    # 	end
    # 	p
    # end;

    # # ╔═╡ bbf648c1-5f6f-456a-960b-4c7374b7c74f
    # p = let
    # 	p = plot(
    # 		p_freq, p_pop_stack;
    # 		layout = grid(1,2),
    # 		size = (1000, 500),
    # 		margin=5mm,
    # 		bottom_margin=20mm,
    # 	)
    # end

    # return p_freq, p_pop_stack
    tvals, Is = compartment_values(solutions)
    Z = sum(Is)
    Is_mut = [sum(Is[findall(idx_mut)[1:n]])./Z for n in 1:N_mut]
    return tvals, Is_mut[end], K

end

# ╔═╡ Cell order:
# ╟─80001173-11b9-4d1c-ab57-ee6edb095fb8
# ╠═a16382a0-8b7b-11ed-050f-b512c0949dd0
# ╟─5094c0ff-bf9e-4476-b0f1-9b9941c1c452
# ╠═446f0274-6944-43a9-bad8-265aec6ada60
# ╠═e1e91203-6a97-43f0-b031-036853f15ee1
# ╠═bbf648c1-5f6f-456a-960b-4c7374b7c74f
# ╠═629df600-203e-4456-b1f0-1479af1a4199
# ╠═3279b428-b03e-48e5-97e8-06f8f8e4d80f
# ╠═30f888a7-d3f4-4089-9c1e-c39ada15444a
# ╠═528bfbbe-0102-447d-a9fd-5e1f17e55238
# ╠═e34bdc45-3bb1-4628-8e55-c40d54d4159a
# ╠═b761969b-8fc1-4772-aeaf-f61c6525fcdc
# ╠═f75d7c6b-ca6f-4dfe-aa00-d460e7776918
# ╠═1237699b-1a27-449d-94c8-983013ce1907
# ╠═0defdd45-4fca-49e4-9216-fbb06d9e0968
# ╟─406c2e67-20e4-4062-aa04-58113b7a06fa
# ╠═9b70216c-13ec-4637-b536-ceb85c7346ab
# ╠═136175fd-5bdf-46ae-8bf4-55d6b9b47f08
# ╠═80fa524d-f0b1-470c-a7e8-f4857e33ee4a
# ╠═424d9f47-c380-4bf4-9dab-be01a4de4830
# ╟─261a3400-1415-4a67-99f7-4d6f2bc2d19d
# ╠═941cb880-c300-4785-bf76-f816aea7bced
# ╠═13b2937c-116e-4476-afe9-5a4fad087c25
# ╠═34c7c53a-dec9-49ae-b2c7-e4a3bcc020c1
# ╟─fc603500-6420-4018-8b5f-a670ef128c5a
# ╟─6571f401-fc4e-4d04-8f8a-83f69b579701
# ╟─36ca0fdb-d0d6-45ee-93f1-090f761e3191
# ╠═41cc2b47-cab0-4f4b-a243-f8adc58deb2d
# ╟─dce2093c-1e9d-47db-8cd1-d47e7a5216cf
# ╟─5e4835fb-3447-4ca4-831b-cabbc2db3528
