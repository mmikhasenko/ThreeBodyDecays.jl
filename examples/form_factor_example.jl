"""
Example demonstrating energy-dependent vertex functions with form-factor functions.

This example shows how to:
1. Create form factors that depend on particle masses
2. Create form factors that depend on energy (invariant mass)
3. Use these form factors in decay chains
"""

using ThreeBodyDecays

# Create a three-body system
masses = ThreeBodyMasses(0.139, 0.139, 0.139; m0 = 0.770)  # πππ from ρ(770)
spins = ThreeBodySpins(1, 1, 0; h0 = 1)  # πππ with ρ(770) intermediate
tbs = ThreeBodySystem(masses, spins)

# Example 1: Mass-dependent form factor
# This form factor depends only on the particle masses
mass_dependent_ff = MassDependentFormFactor((masses) -> 
    masses.m0 / (masses.m1 + masses.m2)  # Simple mass ratio
)

# Example 2: Energy-dependent form factor (Blatt-Weisskopf type)
# This form factor depends on both energy and masses
energy_dependent_ff = EnergyDependentFormFactor((energy, masses) -> 
    1.0 / (1.0 + energy / masses.m0^2)  # Blatt-Weisskopf-like form factor
)

# Example 3: More complex energy-dependent form factor
# This could represent a more sophisticated hadronic form factor
complex_ff = EnergyDependentFormFactor((energy, masses) -> begin
    # Example: form factor that depends on energy and mass differences
    mass_diff = abs(masses.m1 - masses.m2)
    energy_scale = masses.m0^2
    return 1.0 / (1.0 + energy / energy_scale + mass_diff^2 / energy_scale)
end)

# Create vertex functions with different form factors
recoupling = RecouplingLS((2, 2))  # l=1, s=1

# Vertex function with mass-dependent form factor
vf_mass = VertexFunction(recoupling, mass_dependent_ff)

# Vertex function with energy-dependent form factor  
vf_energy = VertexFunction(recoupling, energy_dependent_ff)

# Vertex function with complex form factor
vf_complex = VertexFunction(recoupling, complex_ff)

# Test the form factors with some example values
println("Testing form factors:")
println("Masses: m1=$(masses.m1), m2=$(masses.m2), m3=$(masses.m3), m0=$(masses.m0)")

# Test mass-dependent form factor
mass_factor = mass_dependent_ff.f(masses)
println("Mass-dependent form factor: $mass_factor")

# Test energy-dependent form factors
test_energy = 0.5  # GeV²
energy_factor1 = energy_dependent_ff.f(test_energy, masses)
energy_factor2 = complex_ff.f(test_energy, masses)
println("Energy-dependent form factor (simple): $energy_factor1")
println("Energy-dependent form factor (complex): $energy_factor2")

# Example of how to use in a decay chain
# Note: This is a simplified example - full integration would require
# more complex setup with proper lineshape functions

println("\nForm factors can now be used in decay chains:")
println("- MassDependentFormFactor: depends only on particle masses")
println("- EnergyDependentFormFactor: depends on energy and masses")
println("- Both types can be combined with any recoupling scheme")
println("- Form factors are evaluated during amplitude calculation")
