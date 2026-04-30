using Atlans
using Dates


path_nc = "examples/subsoil-model.nc"
path_csv = "examples/parameters.csv"
forcing = (
    deep_subsidence = DeepSubsidence("examples/subsidence.nc"),
    stage_change = StageChange("examples/change.nc"),
)

model = Model(
    HydrostaticGroundwater,
    DrainingAbcIsotache,
    CarbonStore,
    OverConsolidationRatio,
    AdaptiveCellsize(0.25, 0.01),
    ExponentialTimeStepper(1.0, 2),
    path_nc,
    path_csv,
);

additional_times = map(DateTime, ["2020-01-01", "2025-01-01", "2030-01-01", "2035-01-01"])
simulation = Simulation(
    model,
    "examples/output.nc",
    DateTime("2040-01-01"),
    forcing,
    additional_times,
);
run!(simulation);
