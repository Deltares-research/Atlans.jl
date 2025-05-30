# # Examples of groundwater flow
#
# This example demonstrates two concepts:
#
# * Linear head gradients between boundaries
# * One-dimensional vertical groundwater flow between boundaries

using Revise
using Atlans


const Float = Float

ig = Atlans.InterpolatedGroundwater(
    collect(0.0:1.0:5.0),  # z
    fill(1.0, 5),  # Δz
    Int[1, 5],  # boundary
    Float[3.0, 4.0],  # boundary_ϕ
    fill(false, 5),  # dry
    fill(0.0, 5),  # ϕ
    fill(0.0, 5),  # p
)

function interpolate_head!(ig::Atlans.InterpolatedGroundwater)
    nbound = length(ig.boundary)
    z1 = ig.z[ig.boundary[1]]
    z2 = ig.z[ig.boundary[2]]
    ϕ1 = ig.boundary[1]
    ϕ2 = ig.boundary[2]
    Δz = z2 - z1
    Δϕ = ϕ2 - ϕ1
    j = 2
    z2 = ig.z[ig.boundary[j]]
    ϕ2 = ig.boundary[j]
    i = 1
    while i < nbound
        zi = ig.z[i]
        if zi <= z1
            ig.ϕ[i] = ϕ1
            i += 1
        elseif zi <= z2
            w = (zi - z1) / Δz
            ig.ϕ[i] = ϕ1 + w * Δϕ
            i += 1
        elseif j < nbound
            j += 1
            z1 = z2
            ϕ1 = ϕ2
            z2 = ig.z[ig.boundary[j]]
            ϕ2 = ig.boundary[j]
        else
            break
        end
    end
    ig.ϕ[i:end] .= ϕ2
    return
end

interpolate_head!(ig)


using Plots
