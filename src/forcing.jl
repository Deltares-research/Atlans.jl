using CommonDataModel: @select

abstract type Forcing end


struct Surcharge <: Forcing
	lithology::Array{OptionalInt, 3}
	thickness::Array{OptionalFloat, 3}
	reader::Reader
	lookup::Dict
end


struct StageIndexation <: Forcing
	percentile::Int
	factor::Array{OptionalFloat, 2}
	weir_area::Array{OptionalInt, 2}
	change::Array{OptionalFloat, 2}
	reader::Reader
end


struct DeepSubsidence <: Forcing
	subsidence::Array{OptionalFloat, 2}
	reader::Reader
end


struct StageChange <: Forcing
	change::Array{OptionalFloat, 2}
	reader::Reader
end


struct AquiferHead <: Forcing
	head::Array{OptionalFloat, 2}
	reader::Reader
end


mutable struct Temperature <: Forcing # TODO: maak ruimtelijk
	temp::OptionalFloat
	table::DataFrame
end


struct Forcings
	deep_subsidence::Union{DeepSubsidence, Nothing}
	stage_indexation::Union{StageIndexation, Nothing}
	stage_change::Union{StageChange, Nothing}
	aquifer_head::Union{AquiferHead, Nothing}
	temperature::Union{Temperature, Nothing}
	surcharge::Union{Surcharge, Nothing}
end


function Forcings(;
	deep_subsidence = nothing,
	stage_indexation = nothing,
	stage_change = nothing,
	aquifer_head = nothing,
	temperature = nothing,
	surcharge = nothing,
	bbox::Union{Nothing, NTuple{4, <:Real}} = nothing,
)
	if bbox !== nothing
		deep_subsidence = select_subset(deep_subsidence, bbox)
		stage_indexation = select_subset(stage_indexation, bbox)
		stage_change = select_subset(stage_change, bbox)
		aquifer_head = select_subset(aquifer_head, bbox)
		surcharge = select_subset(surcharge, bbox)
		# temperature is not spatial, so no need to select subset
	end
	Forcings(
		deep_subsidence,
		stage_indexation,
		stage_change,
		aquifer_head,
		temperature,
		surcharge,
	)
end

select_subset(f::Nothing, _) = nothing

Base.size(f::Forcing) = size(f.reader)


function select_subset(r::Reader, bbox::Union{Nothing, NTuple{4, <:Real}})
	xmin, ymin, xmax, ymax = bbox
	subset = @select(r.dataset, $xmin <= x <= $xmax && $ymin <= y <= $ymax)
	return Reader(subset, r.params, r.times)
end


function select_subset(f::StageIndexation, bbox::Union{Nothing, NTuple{4, <:Real}})
	reader = select_subset(f.reader, bbox)
	shape = size(reader)
	return StageIndexation(
		f.percentile,
		fill(1.0, shape),
		Array{OptionalInt}(missing, shape),
		Array{OptionalFloat}(missing, shape),
		reader,
	)
end


function select_subset(f::DeepSubsidence, bbox::Union{Nothing, NTuple{4, <:Real}})
	reader = select_subset(f.reader, bbox)
	shape = size(reader)
	return DeepSubsidence(Array{OptionalFloat}(missing, shape), reader)
end


function select_subset(f::StageChange, bbox::Union{Nothing, NTuple{4, <:Real}})
	reader = select_subset(f.reader, bbox)
	shape = size(reader)
	return StageChange(Array{OptionalFloat}(missing, shape), reader)
end


function select_subset(f::AquiferHead, bbox::Union{Nothing, NTuple{4, <:Real}})
	reader = select_subset(f.reader, bbox)
	shape = size(reader)
	return AquiferHead(Array{OptionalFloat}(missing, shape), reader)
end


function select_subset(f::Surcharge, bbox::Union{Nothing, NTuple{4, <:Real}})
	reader = select_subset(f.reader, bbox)
	shape = size(reader)
	return Surcharge(
		Array{Union{Missing, Int64}}(missing, shape),
		Array{Union{Missing, Float64}}(missing, shape),
		reader,
		f.lookup,
	)
end

function StageIndexation(path::String, percentile::Int)
	reader = prepare_reader(path)
	shape = size(reader)
	return StageIndexation(
		percentile,
		fill(1.0, shape),
		Array{OptionalInt}(missing, shape),
		Array{OptionalFloat}(missing, shape),
		reader,
	)
end


function DeepSubsidence(path::String)
	reader = prepare_reader(path)
	shape = size(reader)
	return DeepSubsidence(Array{OptionalFloat}(missing, shape), reader)
end


function StageChange(path::String)
	reader = prepare_reader(path)
	shape = size(reader)
	return StageChange(Array{OptionalFloat}(missing, shape), reader)
end


function AquiferHead(path::String)
	reader = prepare_reader(path)
	shape = size(reader)
	return AquiferHead(Array{OptionalFloat}(missing, shape), reader)
end


function Surcharge(path_nc::String, path_table::String)
	reader = prepare_reader(path_nc)
	table = prepare_lookup_table(path_table)
	shape = size(reader)

	Surcharge(
		Array{Union{Missing, Int64}}(missing, shape),
		Array{Union{Missing, Float64}}(missing, shape),
		reader,
		table,
	)
end


function Temperature(path::String)
	temperature = CSV.read(path, DataFrame)
	return Temperature(missing, temperature)
end


function read_forcing!(sur::Surcharge, time)
	if time in sur.reader.times
		sur.lithology .= ncread4d(sur.reader, :lithology, time)
		sur.thickness .= ncread4d(sur.reader, :thickness, time)
		return true
	end
	return false
end


function read_forcing!(si::StageIndexation, time)
	if time in si.reader.times
		si.weir_area .= ncread(si.reader, :weir_area, time)
		si.factor .= ncread(si.reader, :factor, time)
		return true
	end
	return false
end


function read_forcing!(ds::DeepSubsidence, time)
	if time in ds.reader.times
		ds.subsidence .= ncread(ds.reader, :subsidence, time)
		return true
	end
	return false
end


function read_forcing!(sc::StageChange, time)
	if time in sc.reader.times
		sc.change .= ncread(sc.reader, :stage_change, time)
		return true
	end
	return false
end


function read_forcing!(ah::AquiferHead, time)
	if time in ah.reader.times
		ah.weir_area .= ncread(ah.reader, :aquifer_change, time)
		return true
	end
	return false
end


function read_forcing!(t::Temperature, time)
	time_idx = findfirst(t.table.times .== time)
	if !isnothing(time_idx)
		t.temp = t.table.temperature[time_idx]
		return true
	end
	return false
end


prepare_forcingperiod!(_::Forcing, _::Model) = nothing


function prepare_forcingperiod!(si::StageIndexation, model::Model)
	change_to_negative = -1
	si.change .= 0.0
	weir_areas = si.weir_area
	replace!(weir_areas, missing => typemin(Int64))
	isarea = fill(false, size(weir_areas))

	for area in unique(weir_areas)
		area == typemin(Int64) && continue

		isarea .= (weir_areas .== area)

		if area == -1 # Change is the subsidence per column
			si.change[isarea] .= model.output.subsidence[isarea] .* change_to_negative
		else
			try
				si.change[isarea] .=
					nanpercentile(model.output.subsidence[isarea], si.percentile) *
					change_to_negative
			catch ArgumentError
				continue
			end
		end
	end
	return
end


function get_elevation_shift(ds::DeepSubsidence, _, I)
	subsidence = ds.subsidence[I]
	ismissing(subsidence) && return 0.0
	return subsidence
end


function get_elevation_shift(si::StageIndexation, _, I)
	change = si.change[I]
	(ismissing(change) || change == 0.0) && return 0.0
	return change
end


function get_elevation_shift(si::StageChange, _, I)
	change = si.change[I]
	(ismissing(change) || change == 0.0) && return 0.0
	return change
end


function apply_forcing!(si::StageIndexation, column, I)
	change = si.change[I]
	factor = si.factor[I]
	(ismissing(change) || change == 0.0 || ismissing(factor)) && return

	set_phreatic_difference!(column, change * factor)
	return
end


function apply_forcing!(ds::DeepSubsidence, column, I)
	subsidence = ds.subsidence[I]
	ismissing(subsidence) && return

	set_deep_subsidence!(column, subsidence)
	return
end


function apply_forcing!(si::StageChange, column, I)
	change = si.change[I]
	(ismissing(change) || change == 0.0) && return

	set_phreatic_difference!(column, change)
	return
end


function apply_forcing!(ah::AquiferHead, column, I)
	head = ah.head[I]
	ismissing(head) && return

	set_aquifer(column, head)
	return
end


function apply_forcing!(t::Temperature, column, _)
	oc = column.oxidation
	for i in 1:length(oc.cells)
		cell = oc.cells[i]
		newcell = update_alpha(cell, t.temp)
		oc.cells[i] = newcell
	end
end


function apply_forcing!(sur::Surcharge, column, I)
	surcharge_column = prepare_surcharge_column(sur, column, I)
	set_surcharge!(column, surcharge_column)
end
