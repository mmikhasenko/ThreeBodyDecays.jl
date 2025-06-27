# ThreeBodyDecays Software Organization

## Overview

The ThreeBodyDecays project is a C++ implementation of three-body decay amplitude calculations, derived from the [ThreeBodyDecays.jl](https://github.com/mmikhasenko/ThreeBodyDecays.jl) Julia library. The software is designed for integration with EvtGen and provides a complete framework for modeling three-body particle decays using partial wave analysis.

## Project Structure

```
ThreeBodyDecays/
├── CMakeLists.txt                     # Build configuration
├── README.md                          # Project documentation
├── build/                             # Build artifacts and executables
├── Examples/                          # Usage examples and demonstrations
├── Tests/                             # Google Test suite
└── Core Files:
    ├── ThreeBodyDecays.{hh,cpp}       # Main decay calculation engine
    ├── ThreeBodyAmplitudeModel.{hh,cpp} # High-level amplitude model manager
    ├── OrientationAngles.{hh,cpp}     # Angular calculations and transformations
    ├── ClebschGordan.{hh,cpp}         # Clebsch-Gordan coefficient calculations
    ├── FormFactors.{hh,cpp}           # Lineshape functions and form factors
    └── jacobi.hh                      # Jacobi polynomial implementation
```

## Core Class Architecture

### 1. Main Library Classes

#### `ThreeBodyDecays` (Core Engine)
- **Purpose**: Central calculation engine for three-body decay amplitudes
- **Key Methods**:
  - `amplitude()`: Calculate decay amplitudes for specific helicity configurations
  - `amplitude4d()`: Calculate full 4D amplitude tensors
  - `intensity()`: Calculate decay intensity/probability density
  - `cosθij()`: Angular calculations between particles
- **Dependencies**: Uses mathematical utilities (ClebschGordan, jacobi.hh, FormFactors) for calculations

#### `ThreeBodyAmplitudeModel` (High-Level Interface)
- **Purpose**: Manages collections of decay chains with coefficients
- **Key Features**:
  - Add/remove decay chains with labels and complex coefficients
  - Calculate total amplitude from multiple contributing channels
  - Compute interference terms between different decay modes
  - Component-wise intensity calculations
- **Dependencies**: Built on top of `DecayChain` and `ThreeBodyDecays`

### 2. Physical System Description Classes

#### `ThreeBodySystem`
- **Purpose**: Describes the physical properties of the three-body system
- **Contains**: 
  - Particle masses (`ms`: array of 4 masses)
  - Spin quantum numbers (`two_js`: doubled spin values)
- **Usage**: Foundation class used by all decay calculations

#### `ThreeBodyParities`
- **Purpose**: Manages parity quantum numbers for all particles
- **Contains**: Parity values ('+' or '-') for each particle
- **Features**: Parity multiplication operations

#### `SpinParity`
- **Purpose**: Parses and manages spin-parity combinations (e.g., "3/2+", "1/2-")
- **Features**: Automatic parsing from string notation

#### `DecayChain`
- **Purpose**: Represents a single decay channel/resonance
- **Contains**:
  - Channel identifier (`k`)
  - Intermediate resonance spin (`two_j`)
  - Lineshape function (`Xlineshape`)
  - Helicity coupling functions (`HRk`, `Hij`)
  - Associated three-body system (`tbs`)

### 3. Angular Mathematics Classes

#### Wigner Rotation Hierarchy
```
AbstractWignerRotation (Abstract Base)
├── TrivialWignerRotation          # Identity rotation
└── WignerRotation (Base Implementation)
    ├── Arg0WignerRotation         # Type 0 angular arguments
    ├── Arg2WignerRotation         # Type 2 angular arguments
    └── Arg3WignerRotation         # Type 3 angular arguments
```

**Purpose**: Calculate Wigner rotation matrices for angular momentum coupling
**Key Features**: 
- Different argument types for different physical configurations
- Cosine calculations for rotation angles (`cos_zeta`)

#### Recoupling Classes
- `NoRecoupling`: Direct helicity coupling without recoupling
- `ParityRecoupling`: Helicity coupling with parity considerations
- `RecouplingType` enum: Specifies recoupling method

### 4. Angular Coordinate System (`OrientationAngles` namespace)

**Note**: The `OrientationAngles` module is independent of the core `ThreeBodyDecays` calculation engine and operates as a separate utility for coordinate transformations.

#### `FourVector`
- **Purpose**: 4-momentum vector operations
- **Features**: Boost transformations, rotations, mass calculations

#### `DecayNode`
- **Purpose**: Represents nodes in decay tree structures
- **Types**: Parent, Resonance, Final particles

#### `HelicityTransformation` & `OrientationAngles`
- **Purpose**: Convert between lab frame and helicity frame coordinates
- **Features**: Angular transformations for partial wave analysis

### 5. Lineshape and Form Factor Classes

#### Lineshape Hierarchy
```
Lineshape (Abstract Base)
├── BreitWigner                    # Standard Breit-Wigner resonance
├── BreitWignerExtended           # Extended BW with form factors
├── MultichannelBreitWigner       # Multi-channel resonances
├── MultichannelBreitWignerMulti  # Advanced multi-channel
├── BuggBW                        # Bugg parameterization
└── Flatte                        # Flatte parameterization
```

**Purpose**: Provide resonance lineshape functions
**Features**: 
- Blatt-Weisskopf form factors
- Breakup momentum calculations
- Various resonance parameterizations

### 6. Mathematical Utilities

#### `ClebschGordan` Module
- **Purpose**: Calculate Clebsch-Gordan coefficients for angular momentum coupling
- **Functions**: 
  - `clebschgordan()`: Standard CG coefficients
  - `CG_doublearg()`: Double-argument version for efficiency
  - Factorial calculations with caching
- **Usage**: Used by `ThreeBodyDecays` for angular momentum coupling calculations

#### `jacobi.hh`
- **Purpose**: Jacobi polynomial calculations (from Boost Math)
- **Usage**: Mathematical foundation for Wigner rotation calculations in `ThreeBodyDecays`
- **Dependencies**: Boost Math library

## Class Interdependencies

### Dependency Flow
```
ThreeBodyAmplitudeModel
    ↓ (uses collections of)
DecayChain
    ↓ (contains)
ThreeBodySystem + LineshapeFunction + HelicityFunction
    ↓ (uses)
ThreeBodyDecays (calculation engine)
    ↓ (uses)
├── Wigner Rotations (AbstractWignerRotation hierarchy)
├── ClebschGordan coefficients
├── jacobi.hh (Jacobi polynomials from Boost Math)
└── FormFactors/Lineshapes

OrientationAngles (INDEPENDENT MODULE)
    ↓ (no dependencies on other core modules)
    ├── FourVector operations
    ├── DecayNode structures
    └── HelicityTransformation calculations
```

### Key Relationships

1. **ThreeBodyAmplitudeModel** → **DecayChain**: One-to-many relationship
2. **DecayChain** → **ThreeBodySystem**: Each chain references the same physical system
3. **DecayChain** → **Lineshape classes**: Function pointer relationships
4. **ThreeBodyDecays** → **Mathematical utilities**: Central calculation engine uses:
   - ClebschGordan coefficients for angular momentum coupling
   - jacobi.hh (Boost Math) for Wigner rotation calculations
   - FormFactors/Lineshapes for resonance parameterizations
5. **Wigner rotations** → **jacobi.hh**: Mathematical dependency for Jacobi polynomial calculations
6. **OrientationAngles** → **Independent**: Self-contained module with no dependencies on core calculation engine

## Build System

### CMake Configuration
- **Library Target**: `ThreeBodyDecaysLib` (static library)
- **Example Executables**: `AmplitudeModelExample`
- **Test Suite**: Google Test framework with individual and combined test runners
- **Compiler Requirements**: C++17 standard

### Testing Structure
- Individual test executables for each major component
- Combined test runner (`RunTests`) for complete validation
- Test coverage includes:
  - Unit tests for mathematical functions
  - Integration tests for amplitude calculations
  - JSON model parsing tests
  - Orientation angle calculations

## Usage Patterns

### Typical Workflow
1. **Define Physical System**: Create `ThreeBodySystem` with masses and spins
2. **Create Decay Chains**: Define `DecayChain` objects with resonance parameters
3. **Build Model**: Use `ThreeBodyAmplitudeModel` to combine multiple chains
4. **Calculate Amplitudes**: Call amplitude methods with kinematic variables
5. **Extract Physics**: Compute intensities, interference terms, etc.

### Module Independence
- **Core Calculation Engine**: `ThreeBodyDecays` with mathematical utilities (ClebschGordan, jacobi.hh, FormFactors)
- **Model Management**: `ThreeBodyAmplitudeModel` for multi-chain calculations
- **Coordinate Transformations**: `OrientationAngles` as independent utility module

### Integration Points
- **EvtGen Integration**: Library designed for use as EvtGen decay model
- **JSON Configuration**: Support for external parameter files
- **Extensibility**: Virtual base classes allow custom lineshape implementations

## Key Design Principles

1. **Modularity**: Clear separation between physics concepts and mathematical implementations
2. **Performance**: Efficient algorithms with caching where appropriate (e.g., factorial calculations)
3. **Flexibility**: Function pointer interfaces allow runtime configuration of lineshapes and couplings
4. **Type Safety**: Strong typing with custom type aliases for physical quantities
5. **Extensibility**: Virtual base classes and template-friendly design for customization

This organization provides a robust foundation for three-body decay modeling while maintaining clear interfaces between different aspects of the calculation.
