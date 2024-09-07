using SafeTestsets
using ThreeBodyDecays

# ThreeBodyDecay structure
@safetestset "Test of the tbs structure" begin
    include("test_three_body_masses.jl")
    include("test_border.jl")
    include("test_dalitz_mapping.jl")
    include("test_three_body_spins.jl")
    include("test_three_body_parities.jl")
    include("test_invariants.jl")
end

@safetestset "Test of the ThreeBodyDecay structure" begin
    include("test_model.jl")
end

@safetestset "Wigner rotation" begin
    include("test_wigner_rotation_dispatch.jl")
    #
    include("test_wigner_permutations.jl")
    include("test_wigner_angle_sumrules.jl")
    #
    include("circle_relations.jl")
    include("test_representation_property.jl")
end

@safetestset "Couplings and jp logic" begin
    include("test_couplings.jl")
    include("test_recoupling.jl")
    include("test_ls_amplitude.jl")
    include("test_sum_over_polarization.jl")
end

@safetestset "Test inference" begin
    include("test_inference.jl")
end

@safetestset "Overall rotation" begin
    include("test_plane_orientation.jl")
end

@safetestset "Integrals" begin
    include("test_integrals.jl")
end
