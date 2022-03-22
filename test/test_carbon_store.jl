@testset "CarbonStore" begin
    f_organic = 0.2
    ρb = 1000.0
    Δz = 1.0
    f_minimum_organic = 0.05
    α = 1.0e-3

    cell = Atlans.CarbonStore(
        Δz,
        f_organic,
        Atlans.mass_organic(f_organic, ρb, Δz),
        Atlans.mass_mineral(f_organic, ρb, Δz),
        Atlans.mass_organic_minimal(
            Atlans.mass_mineral(f_organic, ρb, Δz),
            f_minimum_organic,
        ),
        α,
        0.0,
    )

    @testset "Initialization" begin
        @test typeof(cell) == Atlans.CarbonStore
    end

    @testset "oxidate" begin
        Δt = 1.0
        new_cell = Atlans.oxidate(cell, Δt)
        @test typeof(new_cell) == Atlans.CarbonStore
        @test new_cell.oxidation ≈ 1.0e-7
    end
end
