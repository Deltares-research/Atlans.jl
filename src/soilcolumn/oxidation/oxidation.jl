struct OxidationColumn{O}
    cells::Vector{O}
    z::Vector{Float}
    Δz::Vector{Float}
    result::Vector{Float}
    max_oxidation_depth::Float
    no_oxidation_Δz::Float
    mass_loss::Vector{Float}
end


oxidation_level(column::OxidationColumn{NullOxidation}, _, _, _, _) = nothing
oxidate!(column::OxidationColumn{NullOxidation}, _, _) = nothing
synchronize_z!(column::OxidationColumn{NullOxidation}, _) = nothing


function oxidation_level(
    column::OxidationColumn{O},
    surface_level,
    phreatic_level,
    deep_subsidence,
    phreatic_change,
) where {O<:OxidationProcess}
    oxidation_z = max(
        surface_level - deep_subsidence - column.max_oxidation_depth,
        phreatic_level + phreatic_change + column.no_oxidation_Δz,
    )
    return oxidation_z
end


function oxidate!(
    column::OxidationColumn{O},
    phreatic_level::Float,
    Δt::Float,
) where {O<:OxidationProcess}
    oxidation_z = max( # Doesn't matter if oxidation_z is above surface_level of column
        phreatic_level + column.no_oxidation_Δz,
        surface_level(column) - column.max_oxidation_depth
    )
    column.result .= 0.0
    for index in reverse(1:length(column.cells))
        if column.z[index] > oxidation_z
            cell = column.cells[index]
            newcell = oxidate(cell, Δt)
            column.cells[index] = newcell
            column.result[index] = newcell.oxidation
            column.mass_loss[index] = newcell.Δm
        else
            break
        end
    end
    return
end


function synchronize_z!(column::OxidationColumn{O}, Δz) where {O<:OxidationProcess}
    for i = 1:length(column.cells)
        cell = column.cells[i]
        newcell = @set cell.Δz = Δz[i]
        column.cells[i] = newcell
    end
    return
end
