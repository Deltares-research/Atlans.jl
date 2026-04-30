using Atlans


function carbon_store()
    f_organic = 0.2
    ρb = 1000.0
    Δz = 0.5
    f_minimum_organic = 0.05
    α = 1.0e-3
    return CarbonStore(Δz, f_organic, f_minimum_organic, ρb, α)
end


function consolidation_column(z, Δz)
    cells = fill(NullConsolidation(), length(z))
    σ = fill(NaN, length(z))
    σ′ = fill(NaN, length(z))
    p = fill(NaN, length(z))
    preconsolidation = OverConsolidationRatio(fill(2.15, length(z)))
    result = fill(0.0, length(z))

    return ConsolidationColumn(cells, z, Δz, σ, σ′, p, preconsolidation, result)
end


function oxidation_column(z, Δz)
    cells = fill(NullOxidation(), length(z))
    result = fill(0.0, length(z))

    return OxidationColumn(cells, z, Δz, result, NaN, NaN)
end


function groundwater_column(z)
    phreatic = Phreatic(-1.5)
    dry = fill(false, length(z))
    p = fill(NaN, length(z))

    return HydrostaticGroundwater(z, phreatic, dry, p)
end


function shrinkage_column(z, Δz)
    τ_years = 60.0
    n_vals = [0.7, 0.7, 0.7, 0.7, 0.7, 1.1, 1.2, 1.6, 1.7, 1.7]

    m_clay = 0.8
    m_organic = 0.1

    cells = [SimpleShrinkage(i, n, m_clay, m_organic) for (i, n) in zip(Δz, n_vals)]
    result = fill(NaN, length(z))
    Hv0 = 0.3

    return ShrinkageColumn(cells, z, Δz, result, Hv0)
end


function create_soilcolumn(ncells, thickness, zbase)
    x = 0.0
    y = 0.0

    Δz = fill(thickness, ncells)
    z = (zbase .+ cumsum(Δz)) .- (thickness .* Δz)

    consolidation = consolidation_column(z, Δz) # NullConsolidation
    oxidation = oxidation_column(z, Δz) # NullOxidation
    groundwater = groundwater_column(z)
    shrinkage = shrinkage_column(z, Δz)

    # oxidation = OxidationColumn(
    #     fill(carbon_store(), ncells), z, Δz, fill(0.0, ncells), 1.2
    # )

    return SoilColumn(
        zbase,
        x,
        y,
        z,
        Δz,
        groundwater,
        consolidation,
        oxidation,
        shrinkage,
    )

end


ncells = 10
thickness = 0.5
zbase = -5.0

#%%
ad = AdaptiveCellsize(0.25, 0.01)
timestepper = ExponentialTimeStepper(1.0, 2)
timesteps = create_timesteps(timestepper, 3650.0)

soilcolumn = create_soilcolumn(ncells, thickness, zbase) # soilcolumn with all Atlans attributes

apply_preconsolidation!(soilcolumn)
prepare_forcingperiod!(soilcolumn, 0.01, 0.0, 0.0)
set_phreatic_difference!(soilcolumn, -1.0)

#%%
println(soilcolumn.z)
s, c, o, shr = advance_forcingperiod!(soilcolumn, timesteps)
println(soilcolumn.z)
