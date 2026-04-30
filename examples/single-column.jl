using Revise
using Atlans


function draining_abc_isotache()
    Δz = 0.5
    t = 0.0
    σ′ = NaN
    γ_wet = 12_500.0
    γ_dry = 10_500.0
    c_d = 1.0
    c_v = 0.006912
    U = 0.0
    a = 0.01737
    b = 0.1303
    c = 0.008686
    τ = NaN
    consolidation = NaN

    return DrainingAbcIsotache(
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
    Δz = 0.5
    f_minimum_organic = 0.05
    α = 1.0e-3
    return CarbonStore(Δz, f_organic, f_minimum_organic, ρb, α)
end


function no_oxidation()
    ncell = 20
    thickness = 0.5
    zbase = -10.0
    Δz = fill(thickness, ncell)
    z = zbase .+ cumsum(Δz) .- 0.5 .* Δz

    cc = ConsolidationColumn(
        fill(draining_abc_isotache(), ncell),
        z,
        Δz,
        fill(NaN, ncell), # σ
        fill(NaN, ncell), # σ′
        fill(NaN, ncell), # p
        OverConsolidationRatio(fill(2.15, ncell)),
        fill(NaN, ncell), # result
    )

    oc = OxidationColumn(
        fill(NullOxidation(), ncell),
        z,
        Δz,
        fill(0.0, ncell),
        NaN,
    )

    gw = HydrostaticGroundwater(
        z,
        Phreatic(-0.5),
        fill(false, ncell),
        fill(NaN, ncell),
    )

    timestepper = ExponentialTimeStepper(1.0, 2)
    timesteps = create_timesteps(timestepper, 3650.0)

    column = SoilColumn(zbase, 0.0, 0.0, z, Δz, gw, cc, oc)
    apply_preconsolidation!(column)

    prepare_forcingperiod!(column, 0.0)
    set_phreatic_difference!(column, -1.0)
    for _ in 1:3
        advance_forcingperiod!(column, timesteps)
        prepare_forcingperiod!(column, 0.0)
    end

    @show sum(0.5 .- column.Δz)
    return
end

ncell = 20
thickness = 0.5
zbase = -10.0
Δz = fill(thickness, ncell)
z = zbase .+ cumsum(Δz) .- 0.5 .* Δz

cc = ConsolidationColumn(
    fill(draining_abc_isotache(), ncell),
    z,
    Δz,
    fill(NaN, ncell), # σ
    fill(NaN, ncell), # σ′
    fill(NaN, ncell), # p
    OverConsolidationRatio(fill(2.15, ncell)),
    fill(NaN, ncell), # result
)

oc = OxidationColumn(fill(carbon_store(), ncell), z, Δz, fill(0.0, ncell), 1.2)

gw = HydrostaticGroundwater(
    z,
    Phreatic(-0.5),
    fill(false, ncell),
    fill(NaN, ncell),
)

timestepper = ExponentialTimeStepper(1.0, 2)
timesteps = create_timesteps(timestepper, 3650.0)

column = SoilColumn(zbase, 0.0, 0.0, z, Δz, gw, cc, oc)
apply_preconsolidation!(column)

prepare_forcingperiod!(column, 0.01, 0.0, -1.0)
set_phreatic_difference!(column, -1.0)
s, c, o = advance_forcingperiod!(column, timesteps)
#
for _ in 1:20
    s, c, o = advance_forcingperiod!(column, timesteps)
    prepare_forcingperiod!(column, 0.01, 0.0, 0.0)
end

@show sum(0.5 * 20 - sum(column.Δz))
