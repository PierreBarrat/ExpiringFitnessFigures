### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ 172a8d5a-7404-11ee-1d8f-29a4c743832a
begin
	using Pkg; Pkg.activate("../../../")
	using CSV
	using DataFrames
	using LaTeXStrings
	using Measures
	using StatsBase
	using Plots
end

# ╔═╡ a6beb2f8-afc5-418c-b675-54088e2cd7b2
begin
	dat_low_var = DataFrame(CSV.File("polymorphism_Tc_low_variance/data_panel.csv"))
	transform!(dat_low_var, [:Ne, :av_poly] => ((x1, x2) -> (x2-x1) ./x1) => :rel_err)
	
	dat_high_var = DataFrame(CSV.File("polymorphism_Tc_high_variance/data_panel.csv"))
	transform!(
		dat_high_var, 
		[:Ne, :av_poly] => ((x1, x2) -> (x2-x1) ./x1) => :rel_err
	)
end;

# ╔═╡ 31fba823-c5e4-4090-b32d-57bf9c7dbe6e
let
	p = plot(
		scale=:log10,
		xlabel = L"$Ne = 1/ρ \,\langle β^2 \rangle$",
		ylabel = "Coalescence time",
		frame=:box,
	)
	scatter!(
		dat_low_var.Ne, dat_low_var.av_poly;
		label="Low β variance", marker=(:star, 1, 6),
	)

	scatter!(
		dat_high_var.Ne, dat_high_var.av_poly;
		label="High β variance", marker=(:utriangle, 2, 5),
	)
	
	plot!(
		collect(extrema(dat_low_var.Ne)), collect(extrema(dat_low_var.Ne));
		line=(:black, :dash), label=""
	)

	savefig("Tc_v_Ne.png")
	p
	
end

# ╔═╡ fe599a30-3b4a-411c-97fa-f3b716c1bd45
let
	p = plot(
		# scale=:log10,
		xlabel = "Fraction of sweep overlap",
		ylabel = "\$(T_c - N_e)/N_e\$",
		frame=:box,
	)

	scatter!(
		dat_low_var.frac_overlap, dat_low_var.rel_err;
		label = "Low β variance", marker=(:star, 1, 6)
	)

	scatter!(
		dat_high_var.frac_overlap, dat_high_var.rel_err;
		label = "High β variance", marker=(:utriangle, 1, 5)
	)

	savefig("errorTc_v_sweep_overlap.png")
	p
end

# ╔═╡ 73a12e3a-01dd-462c-b891-fbee3eb9079a
let
	p = plot(
		scale=:log10,
		xlabel = "Ne",
		ylabel = "Coalescence time",
		frame=:box,
		colorbar_title = "Probability of overlap",
		size = (600, 400)
	)
	scatter!(
		dat_low_var.Ne, dat_low_var.av_poly;
		label="Low β variance", marker=:star, markerstrokewidth=0, markersize=7,
		zcolor = dat_low_var.frac_overlap, c=:reds,
	)

	scatter!(
		dat_high_var.Ne, dat_high_var.av_poly;
		label="High β variance", marker=:utriangle, markerstrokewidth=0, markersize=7,
		zcolor = dat_high_var.frac_overlap, c=:reds,
	)
	
	plot!(
		collect(extrema(dat_low_var.Ne)), collect(extrema(dat_low_var.Ne));
		line=(:black, :dash), label=""
	)

	savefig("Tc_v_Ne_colorbar.png")
	p
end

# ╔═╡ Cell order:
# ╠═172a8d5a-7404-11ee-1d8f-29a4c743832a
# ╠═a6beb2f8-afc5-418c-b675-54088e2cd7b2
# ╠═31fba823-c5e4-4090-b32d-57bf9c7dbe6e
# ╟─fe599a30-3b4a-411c-97fa-f3b716c1bd45
# ╠═73a12e3a-01dd-462c-b891-fbee3eb9079a
