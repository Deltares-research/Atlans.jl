const τ_ref = 1.0

abstract type AbcIsotache <: ConsolidationProcess end

struct DrainingAbcIsotache <: AbcIsotache
    Δz::Float64
    t::Float64
    σ′::Float64  # effective stress
    γ_w::Float64  # wet specific mass
    γ_d::Float64  # dry specific mass
    # Degree of consolidation
    c_d::Int64  # drainage one-sided, two-sided
    c_v::Float64  # drainage coefficient
    U::Float64
    # Isotache parameters
    a::Float64
    b::Float64
    c::Float64
    τ::Float64
end

# compute intrinsic time (τ)
function τ_intrinsic(abc::ABC where {ABC<:AbcIsotache}, ocr::Float64)
    if abc.c < 1.0e-4
        return 1.0e-9
    else
        return τ_ref * ocr^((abc.b - abc.a) / abc.c)
    end
end

function τ_intermediate(abc::ABC where {ABC<:AbcIsotache}, loadstep::Float64)
    σ_term = (abc.σ′ - loadstep) / abc.σ′
    abc_term = (abc.b - abc.a) / abc.c
    return abc.τ * σ_term^abc_term
end

function consolidate(
    abc::DrainingAbcIsotache,
    σ′::Float64,
    Δt::Float64,
    )::Tuple{Float64,DrainingAbcIsotache}
    t = abc.t + Δt
    # Degree of consolidation changes
    U = Atlans.U(abc, t) 
    ΔU = U - abc.U 
    # Effective stress changes
    Δσ′ = σ′ - abc.σ′
    σ′ = abc.σ′ + U * Δσ′
    loadstep = ΔU * Δσ′
    # τ changes
    τ_intm = τ_intermediate(abc, loadstep)
    τ = τ_intm + Δt
    # consolidation
    strain = abc.c * log(abc.τ / τ_intm) + log(σ′ / (σ′ - loadstep))
    consolidation = min(abc.Δz, strain * abc.Δz)
    γ_w = Atlans.compress_γ_wet(abc, consolidation)
    γ_d = Atlans.compress_γ_dry(abc, consolidation)
    # return new state
    return consolidation,
    DrainingAbcIsotache(
        abc.Δz - consolidation,  # new
        t,  # new
        σ′, # new
        γ_w,  # new
        γ_d,  # new
        abc.c_d,
        abc.c_v,
        U,  # new
        abc.a,
        abc.b,
        abc.c,
        τ,  # new
    )
end