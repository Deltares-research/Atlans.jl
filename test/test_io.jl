@testset "Input-Output" begin
	using NCDatasets
	using Dates
	using DataFrames

	@testset "read_params_table" begin
		filename = AtlansFixtures.params_table()
		df = Atlans.read_params_table(filename)
		@test typeof(df) == DataFrames.DataFrame
		types = eltype.(eachcol(df))
		@test all(types[1:2] .== String)
		@test all(types[3:4] .== Int)
		@test all(types[5:end] .== Float64)
	end

	@testset "lookup_table" begin
		filename = AtlansFixtures.params_table()
		df = Atlans.read_params_table(filename)
		tables = Atlans.build_lookup_tables(df)
		@test typeof(tables) == Dict{Symbol, Dict{Tuple{Int, Int}, Float64}}
		@test issetequal(
			keys(tables),
			[
				:gamma_wet,
				:gamma_dry,
				:drainage_factor,
				:c_v,
				:a,
				:b,
				:c,
				:ocr,
				:mass_fraction_organic,
				:minimal_mass_fraction_organic,
				:oxidation_rate,
				:rho_bulk,
				:mass_fraction_lutum,
				:shrinkage_degree,
			],
		)
		@test tables[:a][(1, 2)] == 0.01737

		@test Atlans.lookup(tables[:a], [1, 1], [2, 2]) == [0.01737, 0.01737]
	end

	@testset "subsoil_data" begin
		path_csv = AtlansFixtures.params_table()
		path_nc = AtlansFixtures.subsoil_netcdf()
		subsoil = Atlans.prepare_subsoil_data(path_nc, path_csv)

		@test typeof(subsoil) == Atlans.SubsoilData
		@test typeof(subsoil.data) == Dict{Symbol, Array}
		@test Atlans.lookup(subsoil.tables[:a], [1, 1], [2, 2]) == [0.01737, 0.01737]

		lithology = subsoil.data[:lithology]
		@test typeof(lithology) == Array{Int, 3}
		@test size(lithology) == (4, 2, 3)

		phreatic = subsoil.data[:phreatic_level]
		@test typeof(phreatic) == Array{Float64, 2}

		bbox = (30, 30, 70, 70)
		subsoil = Atlans.prepare_subsoil_data(path_nc, path_csv; bbox = bbox)

		@test typeof(subsoil) == Atlans.SubsoilData
		@test typeof(subsoil.data) == Dict{Symbol, Array}
		@test Atlans.lookup(subsoil.tables[:a], [1, 1], [2, 2]) == [0.01737, 0.01737]

		@test issetequal(subsoil.data[:y], [62.5, 37.5])
		@test issetequal(subsoil.data[:x], [37.5])

		lithology = subsoil.data[:lithology]
		@test typeof(lithology) == Array{Int, 3}
		@test size(lithology) == (4, 1, 2)
	end

	@testset "reader" begin
		path_nc = AtlansFixtures.stage_change_netcdf()
		reader = Atlans.prepare_reader(path_nc)

		@test typeof(reader) == Atlans.Reader
		@test typeof(reader.dataset) == NCDatasets.NCDataset{Nothing, Missing}
		@test all(reader.times .== DateTime.(["2020-01-01", "2020-02-01"]))

		diff = Atlans.ncread(reader, :stage_change)
		@test typeof(diff) == Array{Float64, 3}

		diff = Atlans.ncread(reader, :stage_change, DateTime("2020-01-01"))
		@test typeof(diff) == Array{Float64, 2}
	end

	@testset "output_netcdf" begin
		path = tempname()
		x = [12.5, 37.5]
		y = [87.5, 62.5, 37.5]
		ds = Atlans.setup_output_netcdf(path, x, y)

		@test typeof(ds) == NCDatasets.NCDataset{Nothing, Missing}
		@test issetequal(
			keys(ds),
			[
				"time",
				"y",
				"x",
				"phreatic_level",
				"consolidation",
				"oxidation",
				"shrinkage",
				"subsidence",
			],
		)
		@test dimsize(ds["subsidence"]) == (x = 2, y = 3, time = 0)
	end

	@testset "output_writer" begin
		path = tempname()
		x = [12.5, 37.5]
		y = [87.5, 62.5, 37.5]
		writer = Atlans.prepare_writer(path, x, y)

		@test typeof(writer) == Atlans.Writer
		@test typeof(writer.dataset) == NCDatasets.NCDataset{Nothing, Missing}

		index = Atlans.add_time(writer.dataset, DateTime("2020-01-01"))
		@test index == 1

		values = fill(1.0, (2, 3))
		Atlans.ncwrite(writer, :subsidence, values, index)

		index = Atlans.add_time(writer.dataset, DateTime("2020-02-01"))
		@test index == 2

		values = fill(-1.0, (2, 3))
		Atlans.ncwrite(writer, :phreatic_level, values, index)

		@test dimsize(writer.dataset["subsidence"]) == (x = 2, y = 3, time = 2)
		@test dimsize(writer.dataset["phreatic_level"]) == (x = 2, y = 3, time = 2)
	end

	@testset "stage change" begin
		path = AtlansFixtures.stage_change_netcdf()
		forcing = Atlans.StageChange(path)

		@test typeof(forcing) == Atlans.StageChange
		@test all(ismissing.(forcing.change))
		Atlans.read_forcing!(forcing, DateTime("2020-01-01"))
		@test all(forcing.change .≈ -0.1)
	end

	@testset "deep subsidence" begin
		path = AtlansFixtures.deep_subsidence_netcdf()
		forcing = Atlans.DeepSubsidence(path)

		@test typeof(forcing) == Atlans.DeepSubsidence
		@test all(ismissing.(forcing.subsidence))
		Atlans.read_forcing!(forcing, DateTime("2020-01-01"))
		@test all(forcing.subsidence .≈ 0.05)
	end

	@testset "stage indexation" begin
		path = AtlansFixtures.stage_indexation_netcdf()
		forcing = Atlans.StageIndexation(path, 50)

		@test typeof(forcing) == Atlans.StageIndexation
		@test all(ismissing.(forcing.weir_area))
		@test all(forcing.factor .== 1.0)
		Atlans.read_forcing!(forcing, DateTime("2020-01-01"))
		@test all(forcing.weir_area .== 1.0)
		@test all(forcing.factor .== 0.5)
	end

	@testset "surcharge" begin
		path_nc = AtlansFixtures.simple_surcharge_netcdf()
		path_csv = AtlansFixtures.params_table()

		forcing = Atlans.Surcharge(path_nc, path_csv)

		@test typeof(forcing) == Atlans.Surcharge
		@test all(ismissing.(forcing.lithology))
		@test all(ismissing.(forcing.thickness))
		@test size(forcing.lithology) == size(forcing.thickness) == (2, 3, 1)

		Atlans.read_forcing!(forcing, DateTime("2020-01-01"))
		@test all(forcing.lithology .== 2)
		@test all(forcing.thickness .== 0.5)
	end

	@testset "temperature" begin
		path = AtlansFixtures.temperature_table()
		forcing = Atlans.Temperature(path)

		@test typeof(forcing) .== Atlans.Temperature
		@test ismissing(forcing.temp)
		Atlans.read_forcing!(forcing, DateTime("2020-01-01"))
		@test forcing.temp == 14.0
	end

	@testset "forcings" begin
		path_deep_subsidence = AtlansFixtures.deep_subsidence_netcdf()
		path_stage_indexation = AtlansFixtures.stage_indexation_netcdf()
		path_stage_change = AtlansFixtures.stage_change_netcdf()
		path_surcharge = AtlansFixtures.simple_surcharge_netcdf()
		path_temperature = AtlansFixtures.temperature_table()
		path_csv = AtlansFixtures.params_table()

		forcings = Forcings(
			deep_subsidence = DeepSubsidence(path_deep_subsidence),
			stage_indexation = StageIndexation(path_stage_indexation, 50),
			stage_change = StageChange(path_stage_change),
			surcharge = Surcharge(path_surcharge, path_csv),
			temperature = Temperature(path_temperature),
		)
		@test size(forcings.deep_subsidence.reader) == (2, 3)
		@test size(forcings.deep_subsidence) == (2, 3)
		@test size(forcings.stage_indexation) == (2, 3)
		@test size(forcings.stage_change) == (2, 3)
		@test size(forcings.surcharge) == (2, 3, 1)
		@test forcings.aquifer_head === nothing

		# Test select subsets of Forcings
		forcings = Forcings(
			deep_subsidence = DeepSubsidence(path_deep_subsidence),
			stage_indexation = StageIndexation(path_stage_indexation, 50),
			stage_change = StageChange(path_stage_change),
			surcharge = Surcharge(path_surcharge, path_csv),
			temperature = Temperature(path_temperature),
			bbox = (30, 30, 70, 70),
		)
		@test size(forcings.deep_subsidence) == (1, 2)
		@test size(forcings.deep_subsidence.subsidence) == (1, 2)
		@test size(forcings.stage_indexation) == (1, 2)
		@test size(forcings.stage_indexation.factor) == (1, 2)
		@test size(forcings.stage_indexation.weir_area) == (1, 2)
		@test size(forcings.stage_indexation.change) == (1, 2)
		@test size(forcings.stage_change) == (1, 2)
		@test size(forcings.stage_change.change) == (1, 2)
		@test size(forcings.surcharge) == (1, 2, 1)
		@test size(forcings.surcharge.lithology) == (1, 2, 1)
		@test size(forcings.surcharge.thickness) == (1, 2, 1)

		@test all(forcings.deep_subsidence.reader.dataset[:x][:] .== [37.5])
		@test all(forcings.deep_subsidence.reader.dataset[:y][:] .== [62.5, 37.5])
		@test all(forcings.stage_indexation.reader.dataset[:x][:] .== [37.5])
		@test all(forcings.stage_indexation.reader.dataset[:y][:] .== [62.5, 37.5])
		@test all(forcings.stage_change.reader.dataset[:x][:] .== [37.5])
		@test all(forcings.stage_change.reader.dataset[:y][:] .== [62.5, 37.5])
		@test all(forcings.surcharge.reader.dataset[:x][:] .== [37.5])
		@test all(forcings.surcharge.reader.dataset[:y][:] .== [62.5, 37.5])
	end
end
