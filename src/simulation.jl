abstract type Forcing end

struct Reader
    dataset::NCDataset
    params::Dict{Symbol,String}
    times::Vector{DateTime}
end

struct SubsoilData
    data::Dict{Symbol,Array}
    tables::Dict{Symbol,Dict{Tuple{Int,Int},Float}}
end

struct Writer
    dataset::NCDataset
    params::Dict{Symbol,String}
end

struct Output
    x::Vector{Float}
    y::Vector{Float}
    phreatic_level::Array{Float}
    consolidation::Array{Float}
    oxidation::Array{Float}
    shrinkage::Array{Float}
    subsidence::Array{Float}
    carbon_loss::Array{Float}
end


struct Model{G,C,P,O,S,T,A}
    columns::Vector{SoilColumn{G,C,P,O,S}}
    index::Vector{CartesianIndex}
    timestepper::T
    adaptive_cellsize::A
    output::Output
end


struct StageIndexation <: Forcing
    percentile::Int
    factor::Array{OptionalFloat}
    weir_area::Array{OptionalInt}
    change::Array{OptionalFloat}
    reader::Reader
end


struct DeepSubsidence <: Forcing
    subsidence::Array{OptionalFloat}
    reader::Reader
end


struct StageChange <: Forcing
    change::Array{OptionalFloat}
    reader::Reader
end


struct AquiferHead <: Forcing
    head::Array{OptionalFloat}
    reader::Reader
end


mutable struct Temperature <: Forcing # TODO: maak ruimtelijk
    temp::OptionalFloat
    table::DataFrame
end


struct Simulation{F}
    model::Model
    clock::Clock
    writer::Writer
    forcing::F
end

# Initialization
"""
Repeat function not in base / stdlib
"""
function repeat_elements(x, v)
    z = similar(x, sum(v))
    i = 0
    for (x, n) in zip(x, v), _ = 1:n
        @inbounds z[i+=1] = x
    end
    return z
end


"""
In how many equal parts should a thick cell be divided?
"""
function discretize(Δz, maxΔz::Float)
    ncells = convert.(Int, ceil.(Δz ./ maxΔz))
    index = repeat_elements(1:length(Δz), ncells)
    return index, ncells
end


"""
Temporary structure used to create SoilColumns.
"""
struct VerticalDomain
    z::Vector{Float} # cell midpoint
    Δz::Vector{Float}  # cell thickness
    geology::Vector{Int}  # cell geology id
    lithology::Vector{Int}  # cell lithology id
    index::Vector{Int}  # which indices of included column (may contain repeats)
    n::Int  # Total number of cells
end


"""
Temporary structure used to create SoilColumns.
"""
function prepare_domain(
    domainbase,
    modelbase,
    surface,
    thickness,
    Δzmax,
    geology,
    lithology,
)
    is_valid = .!ismissing.(thickness)
    thickness_summed = cumsum(skipmissing(thickness))
    thickness[is_valid] .= thickness_summed

    ztop = modelbase .+ thickness
    zbot = ztop .- thickness

    base_index = findfirst(domainbase .< skipmissing(ztop))
    top_index = findlast(surface .> skipmissing(zbot))

    Δz = thickness[base_index:top_index]
    index, ncells = discretize(Δz, Δzmax)
    Δz = (Δz./ncells)[index]
    z = zbot[base_index] .+ cumsum(Δz) .- 0.5 .* Δz
    n = sum(ncells)
    index .+= (base_index - 1)
    return VerticalDomain(z, Δz, geology[index], lithology[index], index, n)
end


"""
    Model(
        groundwater::Type,
        consolidation::Type,
        oxidation::Type,
        preconsolidation::Type,
        shrinkage::Type,
        adaptive_cellsize,
        timestepper,
        path_subsoil,
        path_lookup
    )

Initialize a model with specified groundwater, consolidation, oxidation and
shrinkage processes from a netCDF file and CSV lookup table describing the
subsoil parameters, appropriate for the chosen processes.
"""
function Model(
    groundwater::Type,
    consolidation::Type,
    oxidation::Type,
    preconsolidation::Type,
    shrinkage::Type,
    adaptive_cellsize,
    timestepper,
    path_subsoil,
    path_lookup,
)
    subsoil = prepare_subsoil_data(path_subsoil, path_lookup)

    x = subsoil.data[:x]
    y = subsoil.data[:y]
    base = subsoil.data[:zbase]
    domainbase = subsoil.data[:domainbase]
    surface = subsoil.data[:surface_level]
    phreatic = subsoil.data[:phreatic_level]
    geology = subsoil.data[:geology]
    lithology = subsoil.data[:lithology]
    thickness = subsoil.data[:thickness]

    columntype = SoilColumn{groundwater,consolidation,preconsolidation,oxidation,shrinkage}
    columns = Vector{columntype}()
    index = Vector{CartesianIndex}()

    for I in CartesianIndices(domainbase)
        (ismissing(domainbase[I]) || ismissing(surface[I]) || ismissing(phreatic[I])) &&
            continue

        domain = prepare_domain(
            domainbase[I],
            base[I],
            surface[I],
            thickness[:, I],
            adaptive_cellsize.Δzmax,
            geology[:, I],
            lithology[:, I],
        )
        length(domain.z) == 0 && continue

        g_column = initialize(groundwater, domain, subsoil, I)
        c_column = initialize(consolidation, preconsolidation, domain, subsoil, I)
        o_column = initialize(oxidation, domain, subsoil, I)
        s_column = initialize(shrinkage, domain, subsoil, I)

        column = SoilColumn(
            domainbase[I],
            x[I[1]],
            y[I[2]],
            domain.z,
            domain.Δz,
            g_column,
            c_column,
            o_column,
            s_column,
        )
        # Set values such as preconsolidation stress, τ0, etc.        
        # This requires groundwater: pore pressure, etc.
        apply_preconsolidation!(column)
        push!(columns, column)
        push!(index, I)
    end

    shape = size(domainbase)
    fillnan() = fill(NaN, shape)

    output = Output(x, y, fillnan(), fillnan(), fillnan(), fillnan(), fillnan(), fillnan())

    return Model(columns, index, timestepper, adaptive_cellsize, output)
end


"""
    Simulation(model, path_output, stop_time)
    
Setup a simulation from an initialized model.
"""
function Simulation(
    model::Model,
    path_output::String,
    stop_time::DateTime,
    forcing,
    additional_times=nothing,
)
    if isnothing(additional_times)
        additional_times = DateTime[]
    end
    writer = prepare_writer(path_output, model.output.x, model.output.y)
    clock = Clock(DateTime[], 1, stop_time)
    simulation = Simulation(model, clock, writer, forcing)
    set_periods!(simulation, additional_times)
    return simulation
end


"""
Advance a single stress period for all columns.

Timesteps are determined by the total duration and the chosen timestepper.
"""
function advance_forcingperiod!(
    model,
    duration;
    deep_subsidence=nothing,
    stage_indexation=nothing,
    stage_change=nothing,
    aquifer_head=nothing,
    temperature=nothing
)
    timesteps = create_timesteps(model.timestepper, duration)
    @progress for (I, column) in zip(model.index, model.columns)
        # Compute pre-loading stresses, set t to 0, etc.
        if isnothing(deep_subsidence)
            column_subsidence = 0.0
        else
            column_subsidence = get_elevation_shift(deep_subsidence, column, I)
        end

        if !isnothing(stage_change)
            column_phreatic_change = get_elevation_shift(stage_change, column, I)
        elseif !isnothing(stage_indexation)
            column_phreatic_change = get_elevation_shift(stage_indexation, column, I)
        else
            column_phreatic_change = 0.0
        end

        prepare_forcingperiod!(
            column,
            model.adaptive_cellsize.split_tolerance,
            column_subsidence,
            column_phreatic_change,
        )
        # Apply changes 
        for forcing in (stage_indexation, deep_subsidence, stage_change, aquifer_head, temperature)
            isnothing(forcing) && continue
            apply_forcing!(forcing, column, I)
        end
        # Compute
        subsidence, consolidation, oxidation, shrinkage, carbon_loss = advance_forcingperiod!(
            column, timesteps
        )
        # Store result
        model.output.phreatic_level[I] = phreatic_level(column.groundwater)

        model.output.subsidence[I] = subsidence
        model.output.consolidation[I] = consolidation
        model.output.oxidation[I] = oxidation
        model.output.shrinkage[I] = shrinkage
        model.output.carbon_loss[I] = carbon_loss
    end
    return
end


function load_forcing!(forcing, key, time, model)
    !haskey(forcing, key) && return nothing
    input = getfield(forcing, key)
    active = read_forcing!(input, time)
    if active
        prepare_forcingperiod!(input, model)
        return input
    else
        return nothing
    end
end


"""
Advance the simulation by a single forcing period.
Reads new input, computes, writes output.
"""
function advance_forcingperiod!(simulation)
    clock = simulation.clock
    time = currenttime(clock)
    duration = periodduration(clock)
    forcing = simulation.forcing
    model = simulation.model

    advance_forcingperiod!(
        model,
        duration;
        deep_subsidence=load_forcing!(forcing, :deep_subsidence, time, model),
        stage_indexation=load_forcing!(forcing, :stage_indexation, time, model),
        stage_change=load_forcing!(forcing, :stage_change, time, model),
        aquifer_head=load_forcing!(forcing, :aquifer_head, time, model),
        temperature=load_forcing!(forcing, :temperature, time, model)
    )

    write(simulation.writer, clock, simulation.model.output)
    advance!(clock)
    return
end


"""
Collect the period boundaries from the forcing input.
"""
function set_periods!(simulation, additional_times)
    clock = simulation.clock
    stop_time = clock.stop_time

    alltimes = DateTime[]
    for forcing in simulation.forcing
        if typeof(forcing) == Temperature
            times = forcing.table.times
        else
            times = forcing.reader.times
        end
        append!(alltimes, times[times.<stop_time])
    end
    for time in additional_times
        time < stop_time && push!(alltimes, time)
    end
    push!(alltimes, stop_time)

    clock.times = sort(unique(alltimes))
    return
end


"""
Run all forcing periods of the simulation.
"""
function run!(simulation)
    clock = simulation.clock

    logger, logfile = init_logger("info", "atlans.log", false)
    with_logger(logger) do
        try
            while currenttime(clock) < clock.stop_time
                advance_forcingperiod!(simulation)
            end
        finally
            close(logfile)
        end
    end
    close(simulation.writer.dataset)
    return
end
