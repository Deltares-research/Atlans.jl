using Test
using Atlans


include("fixtures.jl")

# test cells
include("test_hydrostatic_groundwater.jl")
include("test_carbon_store.jl")
include("test_simple_shrinkage.jl")
#include("test_constant_rate.jl")

# test columns
include("test_draining_abc_isotache.jl")
#include("test_koppejan.jl")
include("test_consolidation_column.jl")
include("test_oxidation_column.jl")
include("test_shrinkage_column.jl")

include("test_timesteppers.jl")
include("test_split.jl")

include("test_column.jl")

include("test_io.jl")
include("test_forcings.jl")
include("test_surcharge.jl")
include("test_simulation.jl")
