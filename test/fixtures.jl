module AtlansFixtures

using Atlans
using Dates
using NCDatasets
using DataFrames
using CSV

function draining_abc_isotache()
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


function carbon_store()
    f_organic = 0.2
    ρb = 1000.0
    Δz = 1.0
    f_minimum_organic = 0.05
    α = 1.0e-3
    return Atlans.CarbonStore(Δz, f_organic, f_minimum_organic, ρb, α)
end

function draining_abc_isotache_column()
    cell = draining_abc_isotache()

    cells = fill(cell, 4)
    z = collect(0.5:1.0:4.0)
    Δz = fill(1.0, 4)
    σ = fill(NaN, 4)
    σ′ = fill(NaN, 4)
    p = fill(NaN, 4)
    result = fill(NaN, 4)
    preconsolidation = Atlans.OverConsolidationRatio(fill(2.15, 4))

    return Atlans.ConsolidationColumn(cells, z, Δz, σ, σ′, p, preconsolidation, result)
end

function carbon_store_column()
    cell = carbon_store()

    cells = fill(cell, 4)
    z = collect(0.5:1.0:4.0)
    Δz = fill(1.0, 4)
    result = fill(NaN, 4)

    return Atlans.OxidationColumn(cells, z, Δz, result)
end

function hydrostatic_groundwater()
    z = collect(0.5:1.0:4.0)
    phreatic = Atlans.Phreatic(3.0)
    dry = fill(false, 4)
    p = fill(NaN, 4)

    return Atlans.HydrostaticGroundwater(z, phreatic, dry, p)
end

function soil_column_hg_abc_cs()
    z = collect(0.5:1.0:4.0)
    Δz = fill(1.0, 4)

    cc = Atlans.ConsolidationColumn(
        fill(draining_abc_isotache(), 4),
        z,
        Δz,
        fill(NaN, 4), # σ
        fill(NaN, 4), # σ′
        fill(NaN, 4), # p
        Atlans.OverConsolidationRatio(fill(2.15, 4)),
        fill(NaN, 4), # result
    )

    oc = Atlans.OxidationColumn(fill(carbon_store(), 4), z, Δz, fill(NaN, 4), 1.0)

    gw =
        Atlans.HydrostaticGroundwater(z, Atlans.Phreatic(3.0), fill(false, 4), fill(NaN, 4))

    return Atlans.SoilColumn(0.0, 0.0, 0.0, z, Δz, gw, cc, oc)
end

function soil_column_hg_abc_null()
    z = collect(0.5:1.0:4.0)
    Δz = fill(1.0, 4)

    cc = Atlans.ConsolidationColumn(
        fill(draining_abc_isotache(), 4),
        z,
        Δz,
        fill(NaN, 4), # σ
        fill(NaN, 4), # σ′
        fill(NaN, 4), # p
        Atlans.OverConsolidationRatio(fill(2.15, 4)),
        fill(NaN, 4), # result
    )

    oc = Atlans.OxidationColumn(fill(Atlans.NullOxidation(), 4), z, Δz, fill(0.0, 4), 1.0)

    gw =
        Atlans.HydrostaticGroundwater(z, Atlans.Phreatic(3.0), fill(false, 4), fill(NaN, 4))

    return Atlans.SoilColumn(0.0, 0.0, 0.0, z, Δz, gw, cc, oc)
end

function soil_column_hg_null_cs()
    z = collect(0.5:1.0:4.0)
    Δz = fill(1.0, 4)

    cc = Atlans.ConsolidationColumn(
        fill(Atlans.NullConsolidation(), 4),
        z,
        Δz,
        fill(NaN, 4), # σ
        fill(NaN, 4), # σ′
        fill(NaN, 4), # p
        Atlans.OverConsolidationRatio(fill(NaN, 4)),
        fill(0.0, 4), # result
    )

    oc = Atlans.OxidationColumn(fill(carbon_store(), 4), z, Δz, fill(NaN, 4), 1.0)

    gw =
        Atlans.HydrostaticGroundwater(z, Atlans.Phreatic(3.0), fill(false, 4), fill(NaN, 4))

    return Atlans.SoilColumn(0.0, 0.0, 0.0, z, Δz, gw, cc, oc)
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
    geology[:] .= 1
    lithology[:] .= 2
    phreatic[:] .= 0.5
    thickness[:] .= 0.25
    base[:] .= 0.0
    domainbase[:] .= 0.0
    surface_level[:] .= 1.0
    max_oxidation_depth .= 1.2
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
        difference[:] .= -0.1
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
        difference[:] .= -0.05
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
    )
    CSV.write(filename, df)
    return filename
end

end # module