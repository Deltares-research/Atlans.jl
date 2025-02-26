module AtlansFixtures

using Atlans
using Dates
using NCDatasets
using DataFrames
using CSV


function draining_abc_isotache_cell()
	Δz = 1.0
	t = 0.0
	σ′ = 10000.0
	γ_wet = 15000.0
	γ_dry = 10000.0
	c_d = 2.0
	c_v = 0.006912
	U = 0.0
	a = 0.01737
	b = 0.1303
	c = 0.008686
	τ = 1.0
	consolidation = 0.0

	return Atlans.DrainingAbcIsotache(
		Δz,
		Δz,
		t,
		σ′,
		γ_wet,
		γ_dry,
		c_d,
		c_v,
		U,
		a,
		b,
		c,
		τ,
		consolidation,
	)
end


function carbon_store_cell()
	f_organic = 0.2
	ρb = 1000.0
	Δz = 1.0
	f_minimum_organic = 0.05
	α = 1.0e-3
	return Atlans.CarbonStore(Δz, f_organic, f_minimum_organic, ρb, α)
end


function shrinkage_cell()
	Δz = 1.0
	n = 1.7
	m_clay = 0.7
	m_organic = 0.2
	return Atlans.SimpleShrinkage(Δz, n, m_clay, m_organic)
end


function draining_abc_isotache_column()
	cell = draining_abc_isotache_cell()

	cells = fill(cell, 4)
	z = collect(0.5:1.0:4.0)
	Δz = fill(1.0, 4)
	σ = fill(NaN, 4)
	σ′ = fill(NaN, 4)
	p = fill(NaN, 4)
	result = fill(NaN, 4)
	preconsolidation = Atlans.OverConsolidationRatio(fill(2.15, 4))

	return Atlans.ConsolidationColumn(
		cells, z, Δz, σ, σ′, p, preconsolidation, result,
	)
end


function carbon_store_column()
	cell = carbon_store_cell()

	cells = fill(cell, 4)
	z = collect(0.5:1.0:4.0)
	Δz = fill(1.0, 4)
	result = fill(NaN, 4)
	carbon_loss = fill(0.0, 4)
	return Atlans.OxidationColumn(cells, z, Δz, result, 1.0, 0.0, carbon_loss)
end


function shrinkage_column()
	cell = shrinkage_cell()

	cells = fill(cell, 4)
	z = collect(0.5:1.0:4.0)
	Δz = fill(1.0, 4)
	result = fill(0.0, 4)

	return Atlans.ShrinkageColumn(cells, z, Δz, result, 0.2)
end


function hydrostatic_groundwater()
	z = collect(0.5:1.0:4.0)
	phreatic = Atlans.Phreatic(3.0)
	dry = fill(false, 4)
	p = fill(NaN, 4)

	return Atlans.HydrostaticGroundwater(z, phreatic, dry, p)
end


"""
Test SoilColumn with all processes: consolidation, oxidation and shrinkage.

"""
function soil_column_hg_abc_cs_shr()
	z = collect(0.5:1.0:4.0)
	Δz = fill(1.0, 4)
	ncells = 4

	cc = Atlans.ConsolidationColumn(
		fill(draining_abc_isotache_cell(), ncells),
		z,
		Δz,
		fill(NaN, ncells), # σ
		fill(NaN, ncells), # σ′
		fill(NaN, ncells), # p
		Atlans.OverConsolidationRatio(fill(2.15, ncells)),
		fill(0.0, ncells), # result
	)

	oc = Atlans.OxidationColumn(
		fill(carbon_store_cell(), ncells),
		z,
		Δz,
		fill(0.0, ncells),
		1.0,
		0.0,
		fill(0.0, ncells),
	)

	shr = Atlans.ShrinkageColumn(
		fill(shrinkage_cell(), ncells), z, Δz, fill(0.0, ncells), 0.2,
	)

	gw = Atlans.HydrostaticGroundwater(
		z, Atlans.Phreatic(3.0), fill(false, ncells), fill(NaN, ncells),
	)

	return Atlans.SoilColumn(0.0, 0.0, 0.0, z, Δz, gw, cc, oc, shr)
end


"""
Test SoilColumn with only consolidation process. Uses NullOxidation and NullShrinkage.

"""
function soil_column_hg_abc_nullox_nullshr()
	z = collect(0.5:1.0:4.0)
	Δz = fill(1.0, 4)
	ncells = 4

	cc = Atlans.ConsolidationColumn(
		fill(draining_abc_isotache_cell(), ncells),
		z,
		Δz,
		fill(NaN, ncells), # σ
		fill(NaN, ncells), # σ′
		fill(NaN, ncells), # p
		Atlans.OverConsolidationRatio(fill(2.15, ncells)),
		fill(0.0, ncells), # result
	)

	oc = Atlans.OxidationColumn(
		fill(Atlans.NullOxidation(), ncells),
		z,
		Δz,
		fill(0.0, ncells),
		1.0,
		0.0,
		fill(0.0, ncells),
	)

	shr = Atlans.ShrinkageColumn(
		fill(Atlans.NullShrinkage(), ncells), z, Δz, fill(0.0, ncells), 0.2,
	)

	gw = Atlans.HydrostaticGroundwater(
		z, Atlans.Phreatic(3.0), fill(false, ncells), fill(NaN, ncells),
	)

	return Atlans.SoilColumn(0.0, 0.0, 0.0, z, Δz, gw, cc, oc, shr)
end


"""
Test SoilColumn with only oxidation process. Uses NullConsolidation and NullShrinkage.

"""
function soil_column_hg_nullcons_cs_nullshr()
	z = collect(0.5:1.0:4.0)
	Δz = fill(1.0, 4)
	ncells = 4

	cc = Atlans.ConsolidationColumn(
		fill(Atlans.NullConsolidation(), ncells),
		z,
		Δz,
		fill(NaN, ncells), # σ
		fill(NaN, ncells), # σ′
		fill(NaN, ncells), # p
		Atlans.OverConsolidationRatio(fill(NaN, ncells)),
		fill(0.0, ncells), # result
	)

	oc = Atlans.OxidationColumn(
		fill(carbon_store_cell(), ncells),
		z,
		Δz,
		fill(0.0, ncells),
		1.0,
		0.0,
		fill(0.0, ncells),
	)

	shr = Atlans.ShrinkageColumn(
		fill(Atlans.NullShrinkage(), ncells), z, Δz, fill(0.0, ncells), 0.2,
	)

	gw = Atlans.HydrostaticGroundwater(
		z, Atlans.Phreatic(3.0), fill(false, ncells), fill(NaN, ncells),
	)

	return Atlans.SoilColumn(0.0, 0.0, 0.0, z, Δz, gw, cc, oc, shr)
end


"""
Test SoilColumn with only shrinkage process. Uses NullConsolidation and NullOxidation.

"""
function soil_column_hg_nullcons_nullox_shr()
	z = collect(0.5:1.0:4.0)
	Δz = fill(1.0, 4)
	ncells = 4

	cc = Atlans.ConsolidationColumn(
		fill(Atlans.NullConsolidation(), ncells),
		z,
		Δz,
		fill(NaN, ncells), # σ
		fill(NaN, ncells), # σ′
		fill(NaN, ncells), # p
		Atlans.OverConsolidationRatio(fill(NaN, ncells)),
		fill(0.0, ncells), # result
	)

	oc = Atlans.OxidationColumn(
		fill(Atlans.NullOxidation(), ncells),
		z,
		Δz,
		fill(0.0, ncells),
		1.0,
		0.0,
		fill(0.0, ncells),
	)

	shr = Atlans.ShrinkageColumn(
		fill(shrinkage_cell(), ncells), z, Δz, fill(0.0, ncells), 0.2,
	)

	gw = Atlans.HydrostaticGroundwater(
		z, Atlans.Phreatic(3.0), fill(false, ncells), fill(NaN, ncells),
	)

	return Atlans.SoilColumn(0.0, 0.0, 0.0, z, Δz, gw, cc, oc, shr)
end


function create_xcoord!(ds, x)
	defVar(
		ds,
		"x",
		x,
		("x",),
		attrib = ["standard_name" => "projection_x_coordinate", "axis" => "X"],
	)
end


function create_ycoord!(ds, y)
	defVar(
		ds,
		"y",
		y,
		("y",),
		attrib = ["standard_name" => "projection_y_coordinate", "axis" => "Y"],
	)
end


function subsoil_netcdf()
	filename = tempname()
	ds = NCDatasets.Dataset(filename, "c")
	defDim(ds, "x", 2)
	defDim(ds, "y", 3)
	defDim(ds, "layer", 4)
	create_xcoord!(ds, [12.5, 37.5])
	create_ycoord!(ds, [87.5, 62.5, 37.5])

	geology = defVar(ds, "geology", Int, ("layer", "x", "y"))
	lithology = defVar(ds, "lithology", Int, ("layer", "x", "y"))
	thickness = defVar(ds, "thickness", Float64, ("layer", "x", "y"))
	phreatic = defVar(ds, "phreatic_level", Float64, ("x", "y"))
	base = defVar(ds, "zbase", Float64, ("x", "y"))
	domainbase = defVar(ds, "domainbase", Float64, ("x", "y"))
	surface_level = defVar(ds, "surface_level", Float64, ("x", "y"))
	max_oxidation_depth = defVar(ds, "max_oxidation_depth", Float64, ("x", "y"))
	hv0_ox = defVar(ds, "no_oxidation_thickness", Float64, ("x", "y"))
	hv0_shr = defVar(ds, "no_shrinkage_thickness", Float64, ("x", "y"))

	geology .= 1
	lithology .= 2
	phreatic .= 0.5
	thickness .= 0.25
	base .= 0.0
	domainbase .= 0.0
	surface_level .= 1.0
	max_oxidation_depth .= 1.2
	hv0_ox .= 0.0
	hv0_shr .= 0.2
	return filename
end


function stage_change_netcdf()
	filename = tempname()
	ds = NCDatasets.Dataset(filename, "c") do ds
		defDim(ds, "x", 2)
		defDim(ds, "y", 3)
		defDim(ds, "time", 2)
		create_xcoord!(ds, [12.5, 37.5])
		create_ycoord!(ds, [87.5, 62.5, 37.5])
		defVar(ds, "time", DateTime.(["2020-01-01", "2020-02-01"]), ("time",))
		difference = defVar(ds, "stage_change", Float64, ("x", "y", "time"))
		difference .= -0.1
	end
	return filename
end


function stage_indexation_netcdf()
	filename = tempname()
	ds = NCDatasets.Dataset(filename, "c") do ds
		defDim(ds, "x", 2)
		defDim(ds, "y", 3)
		defDim(ds, "time", 2)
		create_xcoord!(ds, [12.5, 37.5])
		create_ycoord!(ds, [87.5, 62.5, 37.5])
		defVar(ds, "time", DateTime.(["2020-01-01", "2020-02-01"]), ("time",))
		weir_area = defVar(ds, "weir_area", Float64, ("x", "y", "time"))
		factor = defVar(ds, "factor", Float64, ("x", "y", "time"))

		weir_area .= 1.0
		factor .= 0.5
	end
	return filename
end


function deep_subsidence_netcdf()
	filename = tempname()
	ds = NCDatasets.Dataset(filename, "c") do ds
		defDim(ds, "x", 2)
		defDim(ds, "y", 3)
		defDim(ds, "time", 2)
		create_xcoord!(ds, [12.5, 37.5])
		create_ycoord!(ds, [87.5, 62.5, 37.5])
		defVar(ds, "time", DateTime.(["2020-01-01", "2020-02-01"]), ("time",))
		difference = defVar(ds, "subsidence", Float64, ("x", "y", "time"))
		difference .= 0.05
	end
	return filename
end


function params_table()
	filename = tempname()
	df = DataFrame(
		geology_name = ["NAWA"],
		lithology_name = ["sand"],
		geology = [1],
		lithology = [2],
		gamma_wet = [15000.0],
		gamma_dry = [10000.0],
		drainage_factor = [2.0],
		c_v = [0.006912],
		a = [0.01737],
		b = [0.1303],
		c = [0.008686],
		ocr = [2.15],
		mass_fraction_organic = [0.2],
		minimal_mass_fraction_organic = [0.05],
		oxidation_rate = [0.001],
		rho_bulk = [1000.0],
		mass_fraction_lutum = [0.7],
		shrinkage_degree = [1.4],
	)
	CSV.write(filename, df)
	return filename
end


function temperature_table()
	filename = tempname()
	df = DataFrame(
		times = DateTime.(["2020-01-01", "2020-02-01"]),
		temperature = [14.0, 15.0],
	)
	CSV.write(filename, df)
	return filename
end


function simple_surcharge_netcdf()
	filename = tempname()
	ds = NCDatasets.Dataset(filename, "c") do ds
		defDim(ds, "x", 2)
		defDim(ds, "y", 3)
		defDim(ds, "layer", 1)
		defDim(ds, "time", 1)
		create_xcoord!(ds, [12.5, 37.5])
		create_ycoord!(ds, [87.5, 62.5, 37.5])
		defVar(ds, "layer", [1], ("layer",))
		defVar(ds, "time", DateTime.(["2020-01-01"]), ("time",))
		lithology = defVar(ds, "lithology", Int64, ("x", "y", "layer", "time"))
		thickness = defVar(ds, "thickness", Float64, ("x", "y", "layer", "time"))

		lithology .= 2
		thickness .= 0.5 # 2 cells will be added
	end
	return filename
end


function surcharge_profiles_netcdf()
	filename = tempname()
	ds = NCDatasets.Dataset(filename, "c") do ds
		defDim(ds, "x", 2)
		defDim(ds, "y", 3)
		defDim(ds, "layer", 2)
		defDim(ds, "time", 1)
		create_xcoord!(ds, [12.5, 37.5])
		create_ycoord!(ds, [87.5, 62.5, 37.5])
		defVar(ds, "layer", [1, 2], ("layer",))
		defVar(ds, "time", DateTime.(["2020-01-01"]), ("time",))
		lithology = defVar(ds, "lithology", Int64, ("x", "y", "layer", "time"))
		thickness = defVar(ds, "thickness", Float64, ("x", "y", "layer", "time"))

		lithology .= 2
		thickness[:, :, 1] .= 0.5
		thickness[:, :, 2] .= 0.25
	end
	return filename
end


function test_forcings()
	si = Atlans.StageIndexation(stage_indexation_netcdf(), 50)
	sc = Atlans.StageChange(stage_change_netcdf())
	ds = Atlans.DeepSubsidence(deep_subsidence_netcdf())
	t = Atlans.Temperature(temperature_table())
	sur = Atlans.Surcharge(
		simple_surcharge_netcdf(),
		params_table(),
	)
	return Atlans.Forcings(
		stage_change = sc,
		stage_indexation = si,
		deep_subsidence = ds,
		temperature = t,
		surcharge = sur,
	)
end


function testing_model()
	path_csv = params_table()
	path_nc = subsoil_netcdf()

	return Atlans.Model(
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
end

end # module
