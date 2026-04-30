using Atlans
using Dates

Δzmax = 0.25

basedir = raw"c:/Users/knaake/OneDrive - Stichting Deltares/Documents/bodemdaling/for_deniz"

path_nc = joinpath(basedir, "base_subsurface_model.nc")
path_table = joinpath(basedir, "parameters.csv")

model = Model(
	HydrostaticGroundwater,
	DrainingAbcIsotache,
	CarbonStore,
	OverConsolidationRatio,
	NullShrinkage,
	AdaptiveCellsize(Δzmax, 0.01),
	ExponentialTimeStepper(1.0, 2),
	joinpath(basedir, "base_subsurface_model.nc"),
	joinpath(basedir, "parameters.csv"),
)

forcing = Forcings(
	stage_change = StageChange(joinpath(basedir, "stage_change.nc")),
)

additional_times = map(
	DateTime,
	[
		"1920-01-01",
		"1930-01-01",
		"1940-01-01",
		"1950-01-01",
		"1960-01-01",
		"1970-01-01",
		"1980-01-01",
		"1990-01-01",
		"2000-01-01",
	],
)

simulation = Simulation(
	model,
	joinpath(basedir, "test.nc"),
	DateTime("1940-01-01");
	forcings = forcing,
	additional_times = additional_times,
)

# run!(simulation)
