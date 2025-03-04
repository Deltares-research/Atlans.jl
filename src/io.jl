# Reading
# -------
function prepare_subsoil_data(path_nc, path_csv)
	ds = Dataset(path_nc, "r")
	df = read_params_table(path_csv)
	tables = build_lookup_tables(df)
	data = Dict(Symbol(key) => Array(ds[key]) for key in keys(ds))
	return SubsoilData(data, tables)
end


function prepare_reader(path)
	ds = Dataset(path, "r")
	times = ds["time"][:]
	return Reader(ds, Dict(Symbol(key) => key for key in keys(ds)), times)
end


function ncread(reader, param::Symbol)
	varname = reader.params[param]
	return Array(reader.dataset[varname])
end

function ncread(reader, param::Symbol, t::TimeType)
	varname = reader.params[param]
	times = reader.times
	i = searchsortedfirst(times, t)
	i === nothing && throw(DomainError("time $t not in dataset times $(times)"))
	return reader.dataset[varname][:, :, i]
end


function ncread2d(reader, param::Symbol, index)
	varname = reader.params[param]
	return reader.dataset[varname][index]
end


function ncread3d(reader, param::Symbol, index)
	varname = reader.params[param]
	return reader.dataset[varname][:, index]
end


function ncread4d(reader, param::Symbol, t::TimeType)
	varname = reader.params[param]
	times = reader.times
	i = searchsortedfirst(times, t)
	i === nothing && throw(DomainError("time $t not in dataset times $(times)"))
	return reader.dataset[varname][:, :, :, i]
end


function column_type(_, name)
	if name == :geology || name == :lithology
		return Int
	elseif name == :geology_name || name == :lithology_name
		return String
	else
		return Float
	end
end


function read_params_table(path)
	df = CSV.read(path, DataFrame; delim = ",", stringtype = String, types = column_type)
	return df
end


function build_lookup_tables(df)
	geology = df[:, :geology]
	lithology = df[:, :lithology]
	tables = Dict{Symbol, Dict{Tuple{Int, Int}, Float}}()
	for name in names(df, Float)
		tables[Symbol(name)] = Dict(
			(geo, lit) => val for (geo, lit, val) in zip(geology, lithology, df[:, name])
		)
	end
	return tables
end


function lookup(table, geology, lithology)
	n = length(geology)
	found = Vector{Float}(undef, n)
	for i ∈ 1:n
		found[i] = table[(geology[i], lithology[i])]
	end
	return found
end


function fetch_field(subsoil, field, I, domain)
	if haskey(subsoil.data, field)
		values = subsoil.data[field][domain.index, I]
	elseif haskey(subsoil.tables, field)
		values = lookup(subsoil.tables[field], domain.geology, domain.lithology)
	else
		error("Mandatory field not in netCDF or table: $field")
	end
	return values
end


function fetch_field(table, field, domain)
	if haskey(table, field)
		values = lookup(table[field], domain.geology, domain.lithology)
	else
		error("Mandatory field not in Surcharge lookup table: $field")
	end
	return values
end


fetch_field(reader, ::Type{PreOverburdenPressure}, I, domain) =
	fetch_field(reader, :pop, I, domain)

fetch_field(table, ::Type{PreOverburdenPressure}, domain) =
	fetch_field(table, :pop, domain)

fetch_field(reader, ::Type{OverConsolidationRatio}, I, domain) =
	fetch_field(reader, :ocr, I, domain)

fetch_field(table, ::Type{OverConsolidationRatio}, domain) =
	fetch_field(table, :ocr, domain)


function xy_size(reader)
	size_x = length(reader.dataset["x"])
	size_y = length(reader.dataset["y"])
	return (size_x, size_y)
end


function xyz_size(reader)
	size_x = length(reader.dataset["x"])
	size_y = length(reader.dataset["y"])
	size_z = length(reader.dataset["layer"])
	return (size_x, size_y, size_z)
end

# Writing
# -------
function setup_output_netcdf(path, x, y)::NCDataset
	time_units = "days since 1900-01-01 00:00:00"
	calendar = "proleptic_gregorian"

	ds = NCDatasets.NCDataset(path, "c", format = :netcdf3_64bit_offset)
	defDim(ds, "time", Inf)  # unlimited dimension
	defVar(
		ds,
		"time",
		Float,
		("time",),
		attrib = ["units" => time_units, "calendar" => calendar],
	)
	defVar(
		ds,
		"x",
		x,
		("x",),
		attrib = ["standard_name" => "projection_x_coordinate", "axis" => "X"],
	)
	defVar(
		ds,
		"y",
		y,
		("y",),
		attrib = ["standard_name" => "projection_y_coordinate", "axis" => "Y"],
	)

	defVar(ds, "phreatic_level", Float, ("x", "y", "time"))
	defVar(ds, "consolidation", Float, ("x", "y", "time"))
	defVar(ds, "oxidation", Float, ("x", "y", "time"))
	defVar(ds, "shrinkage", Float, ("x", "y", "time"))
	defVar(ds, "subsidence", Float, ("x", "y", "time"))
	defVar(ds, "carbon_loss", Float, ("x", "y", "time"))
	return ds
end


"Add a new time to the unlimited time dimension, and return the index"
function add_time(ds, time)
	i = length(ds["time"]) + 1
	ds["time"][i] = time
	return i
end


function prepare_writer(path, x, y)
	ds = setup_output_netcdf(path, x, y)
	return Writer(
		ds,
		Dict(
			:phreatic_level => "phreatic_level",
			:consolidation => "consolidation",
			:oxidation => "oxidation",
			:subsidence => "subsidence",
			:shrinkage => "shrinkage",
			:carbon_loss => "carbon_loss",
		),
	)
end


function ncwrite(writer::Writer, param, values, time_index)
	varname = writer.params[param]
	writer.dataset[varname][:, :, time_index] .= values
	return
end


function write(writer, clock, output)
	add_time(writer.dataset, currenttime(clock))
	ncwrite(writer, :phreatic_level, output.phreatic_level, clock.iteration)
	ncwrite(writer, :subsidence, output.subsidence, clock.iteration)
	ncwrite(writer, :consolidation, output.consolidation, clock.iteration)
	ncwrite(writer, :oxidation, output.oxidation, clock.iteration)
	ncwrite(writer, :shrinkage, output.shrinkage, clock.iteration)
	ncwrite(writer, :carbon_loss, output.carbon_loss, clock.iteration)
	return
end
