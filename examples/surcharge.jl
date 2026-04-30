using Atlans
using Dates


function update_domainbase!(subsoil::Atlans.SubsoilData)
    surface = subsoil.data[:surface_level]
    domainbase = subsoil.data[:domainbase]

    pleistocene_to_include = 15 # meter
    domainbase_below_surface = surface .- pleistocene_to_include
    domainbase .-= pleistocene_to_include
    domainbase[ismissing.(domainbase)] = domainbase_below_surface[ismissing.(domainbase)]
end


workdir = raw"c:\Users\knaake\OneDrive - Stichting Deltares\Documents\data\bodemdaling"
path_subsurface_model = joinpath(workdir, "model_gouda.nc")
path_table = joinpath(workdir, "parameters_afwegingskader.csv")
path_surcharge = joinpath(workdir, "4.0m_surcharge.nc")


subsoil = Atlans.prepare_subsoil_data(path_subsurface_model, path_table)

update_domainbase!(subsoil)

model = Model(
    HydrostaticGroundwater,
    DrainingAbcIsotache,
    NullOxidation,
    OverConsolidationRatio,
    NullShrinkage,
    AdaptiveCellsize(0.25, 0.01),
    ExponentialTimeStepper(1.0, 2),
    subsoil,
)

forcing = Forcings(; surcharge = Surcharge(path_surcharge, path_table))

simulation = Simulation(
    model,
    joinpath(workdir, "testrun.nc"),
    DateTime("2080-01-01");
    forcings = forcing,
    additional_times = [DateTime("2023-01-01")],
)

run!(simulation)
