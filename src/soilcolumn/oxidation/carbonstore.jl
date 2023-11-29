using SpecialFunctions: erf


struct CarbonStore <: OxidationProcess
    Δz::Float  # cell thickness
    f_organic::Float  # mass fraction organic material
    f_minimum_organic::Float  # no breakdown beyond this level
    m_organic::Float  # mass of organic material
    m_mineral::Float  # mass of mineral material
    α0::Float # oxidation rate t=0
    α::Float  # oxidation rate
    oxidation::Float  # computed oxidation
    carbon_loss::Float # computed carbon loss of fossil carbon material
end


function CarbonStore(Δz, f_organic, f_minimum_organic, ρb, α)
    m_organic = mass_organic(f_organic, ρb, Δz)
    m_mineral = mass_mineral(f_organic, ρb, Δz)
    return CarbonStore(
        Δz, f_organic, f_minimum_organic, m_organic, m_mineral, α, α, NaN, 0.0
    )
end


function mass_organic(f_organic, ρb, Δz)
    return f_organic * ρb * Δz
end


function mass_mineral(f_organic, ρb, Δz)
    return (1.0 - f_organic) * ρb * Δz
end


function ρ_bulk(m_organic, m_mineral, Δz)
    return (m_organic + m_mineral) / Δz
end


function fraction_organic(m_organic, m_mineral)
    return m_organic / (m_organic + m_mineral)
end


"""
Empirical equation to compute specific volume of organic material.

As the organic matter in a soil breaks down, density increases. The "airiest"
parts are the first to go.
"""
function volume_organic(f_organic, ρb)
    if f_organic == 0
        return 0
    else
        return 0.5 / (f_organic * ρb) * (1.0 + erf((f_organic - 0.2) / 0.1))
    end
end


function fossil_carbon_loss(oxidation, ρb, f_organic)
    fossil_carbon_fraction = 0.55
    oxidation * ρb * f_organic * fossil_carbon_fraction
end


function oxidate(cs::CarbonStore, Δt::Float)
    if cs.α == 0 || cs.Δz == 0.0
        new = @set cs.oxidation = 0.0
        return new
    end
    # Compute new bulk density: may be changed by consolidation
    ρb = ρ_bulk(cs.m_organic, cs.m_mineral, cs.Δz)
    m_minimum_organic = mass_organic(cs.f_minimum_organic, ρb, cs.Δz)
    # Compute change of mass
    Δm = min(cs.m_organic - m_minimum_organic, cs.α * cs.Δz * Δt)
    Δm = max(0, Δm)
    m_organic = cs.m_organic - Δm
    # Convert mass to volume for Δz
    oxidation = volume_organic(cs.f_organic, ρb) * Δm
    f_organic = fraction_organic(m_organic, cs.m_mineral)
    return CarbonStore(
        cs.Δz,
        f_organic,  # new
        cs.f_minimum_organic,
        m_organic, # new
        cs.m_mineral,
        cs.α0,
        cs.α,
        oxidation,  # new
        fossil_carbon_loss(oxidation, ρb, cs.f_organic)
    )
end


"""
    initialize(::Type{CarbonStore}, domain, subsoil, I)

Initialize a OxidationColumn for a domain at location I based subsurface input.
"""
function initialize(::Type{CarbonStore}, domain, subsoil, I)
    f_organic = fetch_field(subsoil, :mass_fraction_organic, I, domain)
    f_minimum_organic = fetch_field(subsoil, :minimal_mass_fraction_organic, I, domain)
    α = fetch_field(subsoil, :oxidation_rate, I, domain)
    ρb = fetch_field(subsoil, :rho_bulk, I, domain)

    max_oxidation_depth = subsoil.data[:max_oxidation_depth][I]
    no_oxidation_Δz = subsoil.data[:no_oxidation_thickness][I]

    cells = Vector{CarbonStore}()
    for (i, Δz) in enumerate(domain.Δz)
        cell = CarbonStore(Δz, f_organic[i], f_minimum_organic[i], ρb[i], α[i])
        push!(cells, cell)
    end

    column = OxidationColumn(
        cells,
        domain.z,
        domain.Δz,
        fill(0.0, domain.n),
        max_oxidation_depth,
        no_oxidation_Δz,
        fill(0.0, domain.n)
    )
    return column
end


"""
    initialize(::Type{SimpleShrinkage}, domain, subsoil, I)

Initialize an empty OxidationColumn (i.e. oxidation is ignored) at location I.
"""
function initialize(::Type{NullOxidation}, domain, _, _)
    cells = Vector{NullOxidation}
    for i in 1:length(domain.Δz)
        push!(cells, NullOxidation())
    end
    return OxidationColumn(
        cells,
        domain.z,
        domain.Δz,
        fill(0.0, domain.n),
        NaN,
        NaN,
        NaN
    )
end
