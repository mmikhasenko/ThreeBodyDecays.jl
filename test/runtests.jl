using SafeTestsets
using ThreeBodyDecays

# ThreeBodyDecay structure
@safetestset "Test of the tbs structure" begin
    include("test_threebodymasses.jl")
    include("test_border.jl")
    include("test_threebodyspins.jl")
    include("test_threebodyparities.jl")
    include("test_invariants.jl")
end

@safetestset "Test of the ThreeBodyDecay structure" begin
    include("test_model.jl")
end

include("wignerrotationdispatch.jl")

# angular functions
@safetestset "Angular Functions block" begin
    include("circle_relations.jl")
    include("test_representaiton_property.jl")
end

@safetestset "Permutations" begin
    include("test_wigner_permutations.jl")
end

@safetestset "Sum rule" begin
    include("wigner_angle_sumrules.jl")
end

# JPC functions
@safetestset "Couplings" begin
    include("test_couplings.jl")
    include("test_recoupling.jl")
    include("test_ls_amplitude.jl")
    include("test_sum_over_polarization.jl")
end

include("test_inference.jl")
