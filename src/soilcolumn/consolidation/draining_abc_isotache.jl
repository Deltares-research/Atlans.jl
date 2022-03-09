struct DrainingAbcIsotache <: AbstractAbcIsotache
    Δz::Float
    Δz_0::Float
    t::Float
    σ′::Float  # effective stress
    γ_wet::Float  # wet specific mass
    γ_dry::Float  # dry specific mass
    # Degree of consolidation
    c_d::Float  # drainage factor
    c_v::Float  # drainage coefficient
    U::Float
    # Isotache parameters
    a::Float
    b::Float
    c::Float
    τ::Float
    consolidation::Float  # Computed consolidation
end

function DrainingAbcIsotache(Δz, γ_wet, γ_dry, c_d, c_v, a, b, c)
    return DrainingAbcIsotache(
        Δz,
        Δz,
        0.0,
        NaN,
        γ_wet,
        γ_dry,
        c_d,
        c_v,
        NaN,
        a,
        b,
        c,
        NaN,
        0.0,
    )
end

function consolidate(abc::DrainingAbcIsotache, σ′::Float, Δt::Float)
    t = abc.t + Δt

    # Degree of consolidation changes
    Unew = U(abc, t)
    ΔU = Unew - abc.U

    # Effective stress changes
    load = σ′ - abc.σ′
    σ′ = abc.σ′ + Unew * load
    loadstep = ΔU * load

    # τ changes
    # This also catches cases where c == 0 for non-creeping soils
    if abc.c < 1.0e-4
        τ⃰ = abc.τ
        τ = τ⃰ + Δt
    else
        τ⃰ = abc.τ * ((abc.σ′ - loadstep) / abc.σ′)^((abc.b - abc.a) / abc.c)
        τ = τ⃰ + Δt
    end

    # consolidation
    elastoplastic = abc.a * log(σ′ / (σ′ - loadstep))
    creep = abc.c * log(τ / τ⃰)
    strain = elastoplastic + creep

    # Thickness should not go below 0
    consolidation = min(abc.Δz, strain * abc.Δz)
    γ_wet = compress_γ_wet(abc, consolidation)
    γ_dry = compress_γ_dry(abc, consolidation)

    # return new state
    return DrainingAbcIsotache(
        abc.Δz - consolidation,  # new
        abc.Δz_0,
        t,  # new
        abc.σ′,
        γ_wet,  # new
        γ_dry,  # new
        abc.c_d,
        abc.c_v,
        Unew,  # new
        abc.a,
        abc.b,
        abc.c,
        τ,  # new,
        consolidation,  # new
    )
end

"""
Turn a collection of vectors into a collection of DrainingAbcIsotache cells.
"""
function draining_abc_isotache_column(
    Δz,
    γ_wet,
    γ_dry,
    c_d,
    c_v,
    a,
    b,
    c,
)::Vector{DrainingAbcIsotache}
    nlayer = length(Δz)
    consolidation = Vector{DrainingAbcIsotache}(undef, nlayer)
    for i = 1:nlayer
        cell = DrainingAbcIsotache(
            Δz[i],
            Δz_0[i],
            0.0,  # t
            0.0,  # σ′
            γ_wet[i],
            γ_dry[i],
            c_d[i],
            c_v[i],
            0.0,  # U
            a[i],
            b[i],
            c[i],
            0.0,  # τ
            0.0,  # consolidation
        )
        consolidation[i] = cell
    end
    return consolidation
end

function initialize(
    ::Type{DrainingAbcIsotache},
    preconsolidation::Type,
    domain,
    subsoil,
    I,
)::ConsolidationColumn{DrainingAbcIsotache}
    γ_wet = fetch_field(subsoil, :gamma_wet, I, domain)
    γ_dry = fetch_field(subsoil, :gamma_dry, I, domain)
    c_d = fetch_field(subsoil, :drainage_factor, I, domain)
    c_v = fetch_field(subsoil, :c_v, I, domain)
    a = fetch_field(subsoil, :a, I, domain)
    b = fetch_field(subsoil, :b, I, domain)
    c = fetch_field(subsoil, :c, I, domain)
    precon_values = fetch_field(subsoil, preconsolidation, I, domain)
    precon = preconsolidation(precon_values)

    cells = Vector{DrainingAbcIsotache}()
    for (i, Δz) in enumerate(domain.Δz)
        cell = DrainingAbcIsotache(Δz, γ_wet[i], γ_dry[i], c_d[i], c_v[i], a[i], b[i], c[i])
        push!(cells, cell)
    end

    z = domain.z
    Δz = domain.Δz
    σ = similar(z)
    σ′ = similar(z)
    p = similar(z)
    result = similar(z)

    return ConsolidationColumn(cells, z, Δz, σ, σ′, p, precon, result)
end

"""
Reset degree of consolidation and time.
"""
function prepare_forcingperiod!(
    column::ConsolidationColumn{DrainingAbcIsotache,P} where {P<:Preconsolidation},
)
    for (i, cell) in enumerate(column.cells)
        cell = @set cell.t = 0.0
        cell = @set cell.U = 0.0
        cell = @set cell.Δz_0 = cell.Δz
        column.cells[i] = cell
    end
    return
end
