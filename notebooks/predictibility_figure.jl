### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ f7892a2c-a108-11ee-28f0-4defac98ddd7
begin
	using Pkg; Pkg.activate("../")
	using CSV
	using Chain
	using DataFrames
	using Distributions
	using LaTeXStrings
	using Measures
	using StatsBase
	using StatsPlots
end

# ╔═╡ 544d0652-cdeb-4495-a8a3-afd9e1a84fdd
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

# ╔═╡ ae30fe90-9aaf-4ef6-add2-87ac34307c06
begin
	plt_defaults = pubfig(24)
	Plots.default(; plt_defaults...)
end

# ╔═╡ edcaa2dc-3715-4d0d-a87d-16595b3064fa
md"# Panel"

# ╔═╡ 916855f7-0a6e-4864-ad74-52c3980fa3f3
pwd()

# ╔═╡ 2c5b4f24-f332-48e9-9e24-3559f4bb29e7
md"## Coalescence time plot"

# ╔═╡ b527df5f-7d64-42b8-89c9-87709bd8ae2d
Tc_folder = "intermediate_notebooks/Tc_from_polymorphism_low_high_variance/"

# ╔═╡ fde7f936-9231-4b2d-b065-f115095e8e73
dat_Tc = let
	df_lowvar = DataFrame(CSV.File(
		joinpath(Tc_folder, "polymorphism_Tc_low_variance/data_panel.csv")
	))	
	df_highvar = DataFrame(CSV.File(
		joinpath(Tc_folder, "polymorphism_Tc_high_variance/data_panel.csv")
	))	
	vcat(df_lowvar, df_highvar)
end;

# ╔═╡ ef39ad11-38d8-4343-8db9-a5cfa2ea006e
plot_Tc = let
	p = plot(
		scale=:log10,
		xlabel = "N_e",
		# ylabel = L"T_{MRCA}",
		ylabel = "coalescence time",
		frame=:box,
		colorbar_title = text("Probability of overlap", 24),
		size = (1200, 800),
		right_margin = 20mm,
		# colorbar=false,
		# xticks = ([10, 100, 1000]),
		lim = (8, 5_000),
		
	)
	scatter!(
		dat_Tc.Ne, dat_Tc.av_poly;
		marker = (:utriangle, 12, stroke(0)),
		label = "",
		zcolor = dat_Tc.frac_overlap, c=:reds,
	)

	# scatter!(
	# 	dat_high_var.Ne, dat_high_var.av_poly;
	# 	label="High β variance", marker=:utriangle, markerstrokewidth=0, markersize=7,
	# 	zcolor = dat_high_var.frac_overlap, c=:reds,
	# )
	
	plot!(
		collect(extrema(dat_Tc.Ne)), collect(extrema(dat_Tc.Ne));
		line=(:black, :dash), label=""
	)

	annotate!(11, 3000, ("E", 50))

	p
end

# ╔═╡ 4ca53bb5-823e-4ab0-9bfd-f20ece94a182
md"## Pfixation plots"

# ╔═╡ 69281b36-2616-4116-adf2-3fa9b120f593
S_to_vec(S) = @chain split(S, r",| |\[|\]") filter!(!isempty, _) parse.(Float64, _)

# ╔═╡ d4f39251-eff9-4164-86da-51dc8e64ba4f
data_pfix = let
	df = DataFrame(CSV.File(
		"intermediate_notebooks/predictability/trajectories_expiring_fitness/data_pfix.csv"
	))
	# array saved to CSV look like "[1,2,3]", function below does the parsing
	select!(
		df, Not(:pfix, :f), 
		[:f, :pfix] => ByRow((x1, x2) -> (S_to_vec(x1), S_to_vec(x2))) => [:f, :pfix]
	)
end;

# ╔═╡ e6b3c4c4-2003-4adc-a1b8-70d3ddfbf4ab
function estimate_pfix(x, s0, α)
	P(β) = (1-β)^(α/s0 - 1)
	eps = 1e-4
	Z = quadgk(β -> P(β), eps, 1-eps)[1]

	return x * quadgk(β -> P(β)/Z, eps, x)[1] + quadgk(β -> β*P(β)/Z, x, 1-eps)[1]
end

# ╔═╡ 795db63d-5612-44ea-9fb9-87abbf99ab69
plots_pfix = let
	αvals = data_pfix.α |> sort |> unique
	ρvals = data_pfix.ρ |> sort |> unique
	s0 = data_pfix.s[1]

	pal = palette(:bluesreds, length(αvals))
	map(enumerate(ρvals)) do (k, ρ)
		dat = sort(data_pfix[data_pfix.ρ .== ρ, :], [:α])
		
		p = plot(
			xlabel = "frequency",
			frame=:box,
			title = "ρ/s = $(ρ/s0)",
			legend = (k == 1 ? :bottomright : false),
			ylabel = (k == 1 ? "P fixation" : ""),
		)
		# annotate!(.8, .05, text("ρ/s = $(ρ/s0)", 24))
		for (i, r) in enumerate(eachrow(dat))
			plot!(
				p, r.f, r.pfix;
				label="α/s=$(round(r.α/r.s; sigdigits=2))", 
				marker=(:cross, 12),
				linewidth = 5,
				color = pal[i],
			)
			annotate!(0.03, 0.94, ("BCD"[k], 50))
		end
		plot!([0,1], [0,1], line=(:black, :dash), label="")
	end
end

# ╔═╡ 0ce58257-c14e-4329-8f3a-df122d5d022f
plots_pfix[3]

# ╔═╡ 73e7f564-c10b-4c48-91d2-8d44676042ed
md"## Inertia plots"

# ╔═╡ 0bd31d06-9430-470b-a49c-27b1e75c2f84
data_inertia = let
	df = DataFrame(CSV.File(
		"intermediate_notebooks/predictability/trajectories_expiring_fitness/data_avtraj.csv"
	))
	select!(
		df, Not(:mt, :mf), 
		[:mt, :mf] => ByRow((x1, x2) -> (S_to_vec(x1), S_to_vec(x2))) => [:mt, :mf]
	)
end;

# ╔═╡ 1412b8a4-d92c-4c39-bf7f-a1ecf5f629dc
ρvals = data_inertia.ρ |> sort |> unique

# ╔═╡ 2a10e3e7-6ae5-461d-a8e9-4100b339b6a7


# ╔═╡ 9b49e765-8b51-4b57-b80d-e9532ea6b964
plot_inertia = let
	αvals = data_inertia.α |> sort |> unique
	ρvals = data_inertia.ρ |> sort |> unique
	s0 = data_inertia.s[1]

	ρ = ρvals[1]
	k = 1
	
	dat = sort(data_inertia[data_inertia.ρ .== ρ, :], [:α])
	pal = palette(:bluesreds, length(αvals)+1)
	
	p = plot(
		xlabel = "time",
		ylabel = "frequency",
		frame=:box,
		title = "",
		legend = (k == 1 ? :bottomright : false),
		xlim = (-100, 500),
		ylim = (-0.01, 1.01),
	)
	
	# annotate!(100, .05, text("ρ/s = $(ρ/s0)", 24))
	
	for (i, r) in enumerate(eachrow(dat))
		plot!(
			p, r.mt, r.mf;
			label="α/s=$(round(r.α/r.s; sigdigits=2))",
			linewidth = 5,
			color = pal[i],
		)
	end

	fb = .5
	df = .01
	hline!([fb], label="", line=(:black, :dash, 8, .25))

	annotate!(-80, 0.94, ("A", 50))
	
	p
end

# ╔═╡ bcff8e16-864e-401d-b419-935499704b3d
let
	l = @layout [
		grid(1,2); grid(1,3)
	]
	p = plot(
		plot_inertia, plot_Tc, plots_pfix...; 
		layout=l, size = (2400, 1600), margin = 15mm,
		dpi = 300,
	)

	savefig("../figures/panel_predictability.png")
	savefig(homedir() * "/Documents/BaleLabo/Notes/ExpiringFitness/figures/panel_predictability.png")
	p
end

# ╔═╡ Cell order:
# ╠═f7892a2c-a108-11ee-28f0-4defac98ddd7
# ╠═544d0652-cdeb-4495-a8a3-afd9e1a84fdd
# ╠═ae30fe90-9aaf-4ef6-add2-87ac34307c06
# ╠═edcaa2dc-3715-4d0d-a87d-16595b3064fa
# ╠═bcff8e16-864e-401d-b419-935499704b3d
# ╠═916855f7-0a6e-4864-ad74-52c3980fa3f3
# ╟─2c5b4f24-f332-48e9-9e24-3559f4bb29e7
# ╠═b527df5f-7d64-42b8-89c9-87709bd8ae2d
# ╠═fde7f936-9231-4b2d-b065-f115095e8e73
# ╠═ef39ad11-38d8-4343-8db9-a5cfa2ea006e
# ╟─4ca53bb5-823e-4ab0-9bfd-f20ece94a182
# ╠═69281b36-2616-4116-adf2-3fa9b120f593
# ╠═d4f39251-eff9-4164-86da-51dc8e64ba4f
# ╠═e6b3c4c4-2003-4adc-a1b8-70d3ddfbf4ab
# ╠═795db63d-5612-44ea-9fb9-87abbf99ab69
# ╠═0ce58257-c14e-4329-8f3a-df122d5d022f
# ╟─73e7f564-c10b-4c48-91d2-8d44676042ed
# ╠═0bd31d06-9430-470b-a49c-27b1e75c2f84
# ╠═1412b8a4-d92c-4c39-bf7f-a1ecf5f629dc
# ╠═2a10e3e7-6ae5-461d-a8e9-4100b339b6a7
# ╠═9b49e765-8b51-4b57-b80d-e9532ea6b964
