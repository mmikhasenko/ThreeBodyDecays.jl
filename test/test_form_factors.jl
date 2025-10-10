using ThreeBodyDecays
using Test

@testset "Form Factor Tests" begin
    # Test basic form factor creation
    @testset "Form Factor Creation" begin
        # Test NoFormFactor (default)
        ff_none = NoFormFactor()
        @test ff_none isa NoFormFactor
        
        # Test MassDependentFormFactor
        ff_mass = MassDependentFormFactor((masses) -> masses.m0 / (masses.m1 + masses.m2))
        @test ff_mass isa MassDependentFormFactor
        
        # Test EnergyDependentFormFactor
        ff_energy = EnergyDependentFormFactor((energy, masses) -> 1.0 / (1.0 + energy / masses.m0^2))
        @test ff_energy isa EnergyDependentFormFactor
    end
    
    @testset "Vertex Function with Form Factors" begin
        # Create a simple recoupling
        recoupling = RecouplingLS((2, 2))  # l=1, s=1
        
        # Test vertex function without form factor
        vf_none = VertexFunction(recoupling)
        @test vf_none.ff isa NoFormFactor
        
        # Test vertex function with mass-dependent form factor
        ff_mass = MassDependentFormFactor((masses) -> masses.m0 / (masses.m1 + masses.m2))
        vf_mass = VertexFunction(recoupling, ff_mass)
        @test vf_mass.ff isa MassDependentFormFactor
        
        # Test vertex function with energy-dependent form factor
        ff_energy = EnergyDependentFormFactor((energy, masses) -> 1.0 / (1.0 + energy / masses.m0^2))
        vf_energy = VertexFunction(recoupling, ff_energy)
        @test vf_energy.ff isa EnergyDependentFormFactor
    end
    
    @testset "Form Factor Evaluation" begin
        # Create test masses
        masses = ThreeBodyMasses(1.0, 2.0, 3.0; m0 = 6.0)
        energy = 4.0  # invariant mass squared
        
        # Test mass-dependent form factor
        ff_mass = MassDependentFormFactor((masses) -> masses.m0 / (masses.m1 + masses.m2))
        result_mass = ff_mass.f(masses)
        @test result_mass ≈ 6.0 / (1.0 + 2.0)  # 2.0
        
        # Test energy-dependent form factor
        ff_energy = EnergyDependentFormFactor((energy, masses) -> 1.0 / (1.0 + energy / masses.m0^2))
        result_energy = ff_energy.f(energy, masses)
        @test result_energy ≈ 1.0 / (1.0 + 4.0 / 36.0)  # 1.0 / (1.0 + 1/9) = 9/10
    end
    
    @testset "Decay Chain with Form Factors" begin
        # Create a simple three-body system
        tbs = ThreeBodySystem(
            ThreeBodyMasses(1.0, 2.0, 3.0; m0 = 6.0),
            ThreeBodySpins(1, 1, 0; h0 = 1)
        )
        
        # Create form factors
        ff_mass = MassDependentFormFactor((masses) -> masses.m0 / (masses.m1 + masses.m2))
        ff_energy = EnergyDependentFormFactor((energy, masses) -> 1.0 / (1.0 + energy / masses.m0^2))
        
        # Create vertex functions with form factors
        vf_mass = VertexFunction(RecouplingLS((2, 2)), ff_mass)
        vf_energy = VertexFunction(RecouplingLS((2, 2)), ff_energy)
        
        # Test that we can create decay chains with form factors
        # This is a basic test - full integration would require more complex setup
        @test vf_mass isa VertexFunction
        @test vf_energy isa VertexFunction
    end
end
