@testset "Simulation" begin

	function fill_optional_float(v, n)
		array = Vector{Union{Float64, Missing}}(undef, n)
		fill!(array, v)
		return array
	end

	function fill_optional_int(v, n)
		array = Vector{Union{Int, Missing}}(undef, n)
		fill!(array, v)
		return array
	end

	@testset "repeat" begin
		out = Atlans.repeat_elements([1.0], 5)
		@test typeof(out) == Array{Float64, 1}
		@test length(out) == 5

		out = Atlans.repeat_elements([1.0, 2.0], [2, 3])
		@test all(out .== [1.0, 1.0, 2.0, 2.0, 2.0])
	end

	@testset "discretize" begin
		Δz = fill(0.5, 4)
		Δz[3] = 0.30
		Δz[4] = 0.20
		out, n = Atlans.discretize(Δz, 0.25)
		@test all(n .== [2, 2, 2, 1])
		@test all(out .== [1, 1, 2, 2, 3, 3, 4])
	end

	@testset "vertical domain" begin
		domainbase             = -5.0
		modelbase              = -10.0
		surface                = 5.0
		Δzmax                 = 0.25
		thickness              = fill_optional_float(0.5, 40)  # from -10.0 to +10.0
		geology                = fill_optional_int(1, 40)
		lithology              = fill_optional_int(2, 40)
		thickness[(end-9):end] .= missing
		geology[(end-9):end]   .= missing
		lithology[(end-9):end] .= missing

		domain = Atlans.prepare_domain(
			domainbase,
			modelbase,
			surface,
			thickness,
			Δzmax,
			geology,
			lithology,
		)

		@test typeof(domain) == Atlans.VerticalDomain
		@test all(domain.z .≈ collect(-4.875:0.25:5.0))
		@test all(domain.Δz .≈ 0.25)
		@test all(domain.geology .== 1)
		@test all(domain.lithology .== 2)
		@test domain.index == Atlans.repeat_elements(11:30, fill(2, 20))
		@test domain.n == 40
		@test domain.zbottom == -5.0
	end

	@testset "vertical domain offsetting domainbase" begin
		modelbase = -2.0
		surface = 0.0
		Δzmax = 0.25
		thickness = fill_optional_float(0.5, 4)
		geology = fill_optional_int(1, 4)
		lithology = fill_optional_int(2, 4)

		offsetting_domainbase = -1.1 # Domainbase depth lies in the middle of a voxel
		domain = Atlans.prepare_domain(
			offsetting_domainbase,
			modelbase,
			surface,
			thickness,
			Δzmax,
			geology,
			lithology,
		)

		@test all(domain.z .≈ collect(-1.375:0.25:0.0))
		@test all(domain.Δz .≈ 0.25)
		@test domain.n == 6
		@test domain.zbottom == -1.5 # Make sure zbottom is derived from cells, not domainbase

		matching_domainbase = -1.0 # Domainbase depth lies at the bottom of a voxel
		domain = Atlans.prepare_domain(
			matching_domainbase,
			modelbase,
			surface,
			thickness,
			Δzmax,
			geology,
			lithology,
		)
		@test domain.n == 4
		@test domain.zbottom == -1.0
	end

	@testset "model" begin
		path_csv = AtlansFixtures.params_table()
		path_nc = AtlansFixtures.subsoil_netcdf()
		timestepper = Atlans.ExponentialTimeStepper(1.0, 2)

		model = Atlans.Model(
			Atlans.HydrostaticGroundwater,
			Atlans.DrainingAbcIsotache,
			Atlans.CarbonStore,
			Atlans.OverConsolidationRatio,
			Atlans.SimpleShrinkage,
			Atlans.AdaptiveCellsize(0.25, 0.01),
			Atlans.ExponentialTimeStepper(1.0, 2),
			path_nc,
			path_csv,
		)

		@test typeof(model) == Atlans.Model{
			Atlans.HydrostaticGroundwater,
			Atlans.DrainingAbcIsotache,
			Atlans.OverConsolidationRatio,
			Atlans.CarbonStore,
			Atlans.SimpleShrinkage,
			Atlans.ExponentialTimeStepper{Int},
			Atlans.AdaptiveCellsize,
		}

		@test length(model.columns) == 6
		@test length(model.index) == 6
	end

	@testset "clock" begin
		clock = Atlans.Clock(
			DateTime.(["2020-01-01", "2020-02-01", "2020-03-01"]),
			1,
			DateTime("2020-03-01"),
		)

		@test Atlans.currenttime(clock) == DateTime("2020-01-01")
		@test Atlans.periodduration(clock) ≈ 31.0
		Atlans.advance!(clock)
		@test Atlans.currenttime(clock) == DateTime("2020-02-01")
		@test Atlans.periodduration(clock) ≈ 29.0
	end

	@testset "simulation" begin
		model = AtlansFixtures.testing_model()
		path_deep_subsidence = AtlansFixtures.deep_subsidence_netcdf()
		deep_subs = DeepSubsidence(path_deep_subsidence)
		forcings = Forcings(deep_subsidence = deep_subs)
		simulation = Simulation(
			model,
			tempname(),
			DateTime("2020-03-01");
			forcings = forcings,
		)

		@test simulation.clock.times ==
			  DateTime.(["2020-01-01", "2020-02-01", "2020-03-01"])
		@test simulation.clock.iteration == 1
		@test simulation.clock.stop_time == DateTime("2020-03-01")

		Atlans.advance_forcingperiod!(simulation)
	end

	@testset "simulation run" begin
		path_csv = AtlansFixtures.params_table()
		path_nc = AtlansFixtures.subsoil_netcdf()
		path_deep_subsidence = AtlansFixtures.deep_subsidence_netcdf()
		path_stage_change = AtlansFixtures.stage_change_netcdf()

		model = Model(
			HydrostaticGroundwater,
			DrainingAbcIsotache,
			CarbonStore,
			OverConsolidationRatio,
			SimpleShrinkage,
			AdaptiveCellsize(0.25, 0.01),
			ExponentialTimeStepper(1.0, 2),
			path_nc,
			path_csv,
		)
		forcings = Forcings(
			deep_subsidence = DeepSubsidence(path_deep_subsidence),
			stage_change = StageChange(path_stage_change),
		)
		outfile = tempname()
		simulation = Simulation(
			model,
			outfile,
			DateTime("2020-03-01");
			forcings = forcings,
		)
		run!(simulation)
		rm("atlans.log")

		ds = Dataset(outfile, "r")
		@test all(keys(ds) .== ["time", "x", "y", "phreatic_level", "consolidation", "oxidation", "shrinkage", "subsidence"])
		@test size(ds["subsidence"]) == (2, 3, 2)
	end

	@testset "simulation run subset" begin
		path_csv = AtlansFixtures.params_table()
		path_nc = AtlansFixtures.subsoil_netcdf()
		path_deep_subsidence = AtlansFixtures.deep_subsidence_netcdf()
		path_stage_change = AtlansFixtures.stage_change_netcdf()
		subset_bbox = (30, 30, 70, 70)

		model = Model(
			HydrostaticGroundwater,
			DrainingAbcIsotache,
			CarbonStore,
			OverConsolidationRatio,
			SimpleShrinkage,
			AdaptiveCellsize(0.25, 0.01),
			ExponentialTimeStepper(1.0, 2),
			path_nc,
			path_csv;
			bbox = subset_bbox,
		)
		forcings = Forcings(
			deep_subsidence = DeepSubsidence(path_deep_subsidence),
			stage_change = StageChange(path_stage_change),
			bbox = subset_bbox,
		)
		outfile = tempname()
		simulation = Simulation(
			model,
			outfile,
			DateTime("2020-03-01");
			forcings = forcings,
		)
		run!(simulation)
		rm("atlans.log")

		ds = Dataset(outfile, "r")
		@test all(keys(ds) .== ["time", "x", "y", "phreatic_level", "consolidation", "oxidation", "shrinkage", "subsidence"])
		@test size(ds["subsidence"]) == (1, 2, 2)
	end
end
