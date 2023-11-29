mutable struct Base
    z::Float
end

"""
x, y, z are all midpoints.
"""
struct SoilColumn{G,C,P,O,S}
    base::Base
    x::Float
    y::Float
    z::Vector{Float}
    Δz::Vector{Float}
    groundwater::G
    consolidation::ConsolidationColumn{C,P}
    oxidation::OxidationColumn{O}
    shrinkage::ShrinkageColumn{S}
    subsidence::Vector{Float}
end

function SoilColumn(
    base::Float,
    x,
    y,
    z,
    Δz,
    groundwater,
    consolidation,
    oxidation,
    shrinkage,
)
    return SoilColumn(
        Base(base),
        x,
        y,
        z,
        Δz,
        groundwater,
        consolidation,
        oxidation,
        shrinkage,
        fill(NaN, length(z)),
    )
end

surface_level(column) = column.z[end] + 0.5 * column.Δz[end]

# Forcing

#struct Surcharge{G,C,O}
#    Δz::Vector{Float}
#    groundwater::G
#    consolidation::ConsolidationColumn{C}
#    oxidation::GroundwaterColumn{O}
#end
#
#function set_surcharge!(
#    column::SoilColumn,
#    surcharge,
#)
#    set_surcharge!(column.groundwater, surcharge.groundwater)
#    set_surcharge!(column.consolidation, surcharge.consolidation)
#    set_surcharge!(column.oxidation, surcharge.oxidation)
#    return
#end

function set_deep_subsidence!(column::SoilColumn, subsidence::Float)
    column.base.z -= subsidence
    update_z!(column)
end

function set_aquifer!(column::SoilColumn, ϕ)
    set_aquifer!(column.groundwater, ϕ)
    return
end

function set_aquifer_difference!(column::SoilColumn, Δϕ)
    set_aquifer_difference!(column.groundwater, Δϕ)
    return
end

function set_phreatic!(column::SoilColumn, ϕ)
    set_phreatic!(column.groundwater, ϕ)
    return
end

function set_phreatic_difference!(column::SoilColumn, Δϕ)
    set_phreatic_difference!(column.groundwater, Δϕ)
end

# Computation

function exchange_pore_pressure!(column::SoilColumn)
    column.consolidation.p .= column.groundwater.p * γ_water
end

"""
    prepare_timestep!(column)
    
Prepare a single timestep. This updates stresses in the column and accounts for
e.g. drowning of the column.

This computes:

    * Pore pressure
    * Total stress
    * Effective stress
"""
function prepare_timestep!(column::SoilColumn, Δt)
    flow!(column.groundwater, Δt)
    exchange_pore_pressure!(column)
    total_stress!(column.consolidation, phreatic_level(column.groundwater))
    effective_stress!(column.consolidation)
end

function initial_stress!(column)
    prepare_timestep!(column, 0.0)
    effective_stress!(column.consolidation)
    transfer_stress!(column.consolidation)
    return
end

function apply_preconsolidation!(column)
    initial_stress!(column)
    apply_preconsolidation!(column.consolidation)
end

"""
    prepare_forcingperiod!(column, split_tolerance)
    
Prepare a single forcing period.

This:

    * splits the soil column at maximum oxidation depth (if necessary)
    * computes pore pressure, total stress, effective stress prior to this stress periods loading
    * sets the effective stress in every consolidation cell
    * Reset U and t for DrainingConsolidation processes.
    
Note that splitting requires knowing where the phreatic level ends up after
applying all forcings. This means that changes to phreatic level and deep
subsidence must be accounted for. Furthermore, the split must be applied
before applying the changes to compute the pre-loading effective stress.
"""
function prepare_forcingperiod!(
    column::SoilColumn,
    split_tolerance,
    deep_subsidence=0.0,
    phreatic_change=0.0,
)
    surface = surface_level(column)
    phreatic = phreatic_level(column.groundwater)

    # First split for oxidation process
    oxidation_z = oxidation_level(
        column.oxidation,
        surface,
        phreatic,
        deep_subsidence,
        phreatic_change,
    )
    if !isnothing(oxidation_z)
        split!(column, oxidation_z, split_tolerance)
    end

    # Now split for shrinkage process
    shrinkage_z = shrinkage_level(
        column.shrinkage,
        phreatic,
        phreatic_change,
    )
    if !isnothing(shrinkage_z)
        split!(column, shrinkage_z, split_tolerance)
    end

    initial_stress!(column)
    prepare_forcingperiod!(column.consolidation)
    return
end


"""
    function update_z!(column)
    
Compute new midpoints and surface level.
"""
function update_z!(column::SoilColumn)
    column.z .= column.base.z .+ cumsum(column.Δz) .- 0.5 .* column.Δz
end


"""
    subside!(column)
    
Apply consolidation, oxidation and shrinkage to thickness
"""
function subside!(column::SoilColumn)
    # Δz should not become negative
    column.subsidence .= min.(
        (
            column.consolidation.result
            .+
            column.oxidation.result
            .+
            column.shrinkage.result
        ),
        column.Δz,
    )
    column.Δz .-= column.subsidence
    synchronize_z!(column.groundwater, column.Δz)
    synchronize_z!(column.consolidation, column.Δz)
    synchronize_z!(column.oxidation, column.Δz)
    synchronize_z!(column.shrinkage, column.Δz)
    update_γ!(column.consolidation, column.shrinkage.result)
    update_z!(column)
end

"""
    advance_timestep!(column, Δt)

Advance a single timestep.

During a timestep the following states are computed:

    * head
    * pore pressure
    * total stress
    * effective stress
    
Then, consolidation and oxidation are computed.

Finally, thickness and elevation are updated.
"""
function advance_timestep!(c::SoilColumn, Δt::Float)
    prepare_timestep!(c, Δt)
    consolidate!(c.consolidation, phreatic_level(c.groundwater), Δt)
    oxidate!(c.oxidation, phreatic_level(c.groundwater), Δt)
    shrink!(c.shrinkage, phreatic_level(c.groundwater), Δt)
    subside!(c)
    return (
        sum(c.subsidence),
        sum(c.consolidation.result),
        sum(c.oxidation.result),
        sum(c.shrinkage.result),
        sum(c.oxidation.carbon_loss)
    )
end

"""
    advance_forcingperiod!(column, timesteps)
    
Advances a prepared column by a number of timesteps.

Note, the correct order of execution is:

* prepare a forcing period: compute pre-load stress
* apply forcing: change load
* advance forcing period

"""
function advance_forcingperiod!(c::SoilColumn, timesteps::Vector{Float})
    subsidence = 0.0
    consolidation = 0.0
    oxidation = 0.0
    shrinkage = 0.0
    carbon_loss = 0.0
    for Δt in timesteps
        Δs, Δc, Δo, Δsh, Δm = advance_timestep!(c, Δt)
        subsidence += Δs
        consolidation += Δc
        oxidation += Δo
        shrinkage += Δsh
        carbon_loss += Δm
    end
    return subsidence, consolidation, oxidation, shrinkage, carbon_loss
end

# Output
output(oc::OxidationColumn) = sum(cell.oxidation for cell in oc.cells)
output(cc::ConsolidationColumn) = sum(cell.consolidation for cell in cc.cells)
output(gw::GW where {GW<:GroundwaterColumn}) = phreatic_level(gw)
output(sc::ShrinkageColumn) = sum(cell.shrinkage for cell in sc.cells)
