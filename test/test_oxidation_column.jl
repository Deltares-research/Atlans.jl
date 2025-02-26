@testset "OxidationColumn" begin
    column = AtlansFixtures.carbon_store_column()

    @testset "constructor" begin
        @test typeof(column) == Atlans.OxidationColumn{
            Atlans.CarbonStore
        }
    end

    @testset "surface_level" begin
        @test Atlans.surface_level(column) ≈ 4.0
    end

    @testset "oxidation_depth" begin
        phreatic = 2.5
        oz = max(
            phreatic + column.no_oxidation_Δz,
            Atlans.surface_level(column) - column.max_oxidation_depth
        )
        @test oz ≈ 3.0
    end

    @testset "oxidate" begin
        phreatic = 2.5
        Δt = 3650.0

        Atlans.oxidate!(column, phreatic, Δt)

        expected_ox = [0.0, 0.0, 0.0, 9.125e-3]
        expected_forg = [0.2, 0.2, 0.2, 0.19706930295578862]
        expected_carbon_loss = [0.0, 0.0, 0.0, 1.0037500000000001]

        @test all(column.result .≈ expected_ox)
        @test all(
            c.f_organic ≈ forg for (c, forg) in zip(column.cells, expected_forg)
        )
        @test all(column.carbon_loss .≈ expected_carbon_loss)
    end
end
