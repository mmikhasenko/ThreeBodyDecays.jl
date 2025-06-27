# ThreeBodyDecays Package Dependency Diagram

## Overview
This diagram shows the isolated, self-contained packages and the larger classes that depend on them in the ThreeBodyDecays library.

## Package Architecture

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           ISOLATED PACKAGES                                 │
│                    (Self-contained, no external dependencies)               │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│  📦 ClebschGordan Package                                                  │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │ • clebschgordan() - CG coefficients                                │   │
│  │ • CG_doublearg() - Double-argument version                         │   │
│  │ • getLogFactorialcg() - Factorial calculations                     │   │
│  │ • f_logfact2() - Log factorial utilities                           │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│  📦 jacobi.hh Package                                                     │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │ • boost::math::jacobi<double>() - Jacobi polynomials               │   │
│  │ • Mathematical foundation for Wigner rotations                     │   │
│  │ • External dependency: Boost Math library                          │   │
│  │ • External dependency: Boost Math library                          │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│  📦 FormFactors Package                                                   │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │ • BlattWeisskopf() - Form factor calculations                      │   │
│  │ • breakup() - Breakup momentum calculations                        │   │
│  │ • Lineshape classes:                                               │   │
│  │   - BreitWigner                                                    │   │
│  │   - MultichannelBreitWigner                                        │   │
│  │   - BuggBW                                                         │   │
│  │   - Flatte                                                         │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│  📦 OrientationAngles Package (INDEPENDENT)                              │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │ • FourVector - 4-momentum operations                               │   │
│  │ • DecayNode - Decay tree structures                                │   │
│  │ • HelicityTransformation - Coordinate transformations              │   │
│  │ • OrientationAngles - Angular calculations                         │   │
│  │ • NO DEPENDENCIES on core calculation engine                       │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│                           CORE CLASSES                                     │
│                    (Depend on isolated packages)                           │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│  🔧 ThreeBodyDecays (Core Calculation Engine)                            │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │ DEPENDENCIES:                                                       │   │
│  │ ↓ ClebschGordan Package                                             │   │
│  │ ↓ jacobi.hh Package                                                 │   │
│  │ ↓ FormFactors Package                                               │   │
│  │                                                                     │   │
│  │ KEY METHODS:                                                        │   │
│  │ • amplitude() - Single chain calculations                          │   │
│  │ • amplitude4d() - 4D tensor calculations                           │   │
│  │ • intensity() - Probability density                                │   │
│  │ • cosθij() - Angular calculations                                  │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│  🏗️ ThreeBodyAmplitudeModel (High-Level Manager)                         │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │ DEPENDENCIES:                                                       │   │
│  │ ↓ ThreeBodyDecays (uses internally)                                │   │
│  │ ↓ DecayChain (manages collections)                                 │   │
│  │                                                                     │   │
│  │ KEY METHODS:                                                        │   │
│  │ • add() - Add decay chains with coefficients                       │   │
│  │ • amplitude4d() - Total amplitude from all chains                 │   │
│  │ • intensity() - Total intensity with interference                  │   │
│  │ • interference_terms() - Cross-channel interference                │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│                           SUPPORTING CLASSES                               │
│                    (Used by core classes)                                  │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│  📊 Physical System Classes                                               │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │ • ThreeBodySystem - Masses and spins                               │   │
│  │ • ThreeBodyParities - Parity quantum numbers                       │   │
│  │ • SpinParity - Spin-parity parsing                                 │   │
│  │ • DecayChain - Single decay channel                                │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│  🔄 Wigner Rotation Classes                                               │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │ • AbstractWignerRotation (base)                                     │   │
│  │ • TrivialWignerRotation                                             │   │
│  │ • WignerRotation (base)                                             │   │
│  │ • Arg0WignerRotation                                                │   │
│  │ • Arg2WignerRotation                                                │   │
│  │ • Arg3WignerRotation                                                │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│  🔗 Recoupling Classes                                                    │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │ • NoRecoupling                                                      │   │
│  │ • ParityRecoupling                                                  │   │
│  │ • RecouplingType enum                                               │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────────────────────┘

## Dependency Flow Diagram

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           ISOLATED PACKAGES                                 │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────┐    ┌─────────────┐    ┌─────────────┐    ┌─────────────┐
│ClebschGordan│    │  jacobi.hh  │    │FormFactors  │    │Orientation  │
│   Package   │    │   Package   │    │  Package    │    │  Angles     │
│             │    │             │    │             │    │  Package    │
│ • CG coeffs │    │ • Jacobi    │    │ • Lineshapes│    │ • 4-Vectors │
│ • Factorials│    │ • Wigner    │    │ • Form      │    │ • Decay     │
│ • Utilities │    │   support   │    │   factors   │    │   trees     │
└─────────────┘    └─────────────┘    └─────────────┘    └─────────────┘
       │                   │                   │                   │
       │                   │                   │                   │
       │                   │                   │                   │ (INDEPENDENT)
       │                   │                   │                   │
       ▼                   ▼                   ▼                   │
┌─────────────────────────────────────────────────────────────────────────────┐
│                           CORE CALCULATION ENGINE                           │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│                    ThreeBodyDecays                                         │
│                                                                             │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │ • Uses ClebschGordan for angular momentum coupling                  │   │
│  │ • Uses jacobi.hh for Wigner rotation calculations                   │   │
│  │ • Uses FormFactors for resonance lineshapes                         │   │
│  │ • Contains Wigner rotation hierarchy                                │   │
│  │ • Contains Physical system classes                                  │   │
│  │ • Contains Recoupling classes                                       │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────────────────────┘
                                │
                                ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                        HIGH-LEVEL MODEL MANAGER                             │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│                ThreeBodyAmplitudeModel                                     │
│                                                                             │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │ • Uses ThreeBodyDecays internally                                   │   │
│  │ • Manages collections of DecayChain objects                         │   │
│  │ • Calculates interference between multiple channels                 │   │
│  │ • Provides high-level interface for complete decay models          │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────────────────────┘
```

## Key Architectural Principles

### 1. **Isolation of Mathematical Packages**
- **ClebschGordan**: Pure mathematical functions, no external dependencies
- **jacobi.hh**: Wrapper around Boost Math, minimal interface
- **FormFactors**: Self-contained lineshape implementations
- **OrientationAngles**: Completely independent coordinate transformation utility

### 2. **Layered Architecture**
- **Isolated Packages** → **Core Engine** → **High-Level Manager**
- Each layer depends only on the layer below it
- Clear separation of concerns

### 3. **Dependency Minimization**
- Core calculation engine (`ThreeBodyDecays`) uses isolated packages
- High-level manager (`ThreeBodyAmplitudeModel`) uses core engine
- No circular dependencies
- `OrientationAngles` remains completely independent

### 4. **Extensibility**
- New lineshapes can be added to `FormFactors` package
- New mathematical utilities can be added as isolated packages
- Core engine can be extended without affecting high-level interface

## Usage Patterns

### **Direct Package Usage**
```cpp
// Use isolated packages directly
auto cg = clebschgordan(j1, m1, j2, m2, j, m);
auto jacobi = boost::math::jacobi<double>(n, a, b, z);
auto bw = FormFactors::make_breit_wigner(mass, width);
```

### **Core Engine Usage**
```cpp
// Use core calculation engine
ThreeBodyDecays tbd;
auto amplitude = tbd.amplitude(decay_chain, σs, helicities, k_amp);
```

### **High-Level Model Usage**
```cpp
// Use high-level model manager
ThreeBodyAmplitudeModel model;
model.add(decay_chain1, "rho", 1.0);
model.add(decay_chain2, "f2", 0.5);
auto total_intensity = model.intensity(σs, k_amp);
```

This architecture provides a robust, maintainable, and extensible foundation for three-body decay calculations while maintaining clear boundaries between different levels of abstraction. 