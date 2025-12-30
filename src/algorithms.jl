# ODEInterface.jl Algorithms

abstract type ODEInterfaceAlgorithm <: DiffEqBase.AbstractODEAlgorithm end
abstract type ODEInterfaceImplicitAlgorithm <: ODEInterfaceAlgorithm end
abstract type ODEInterfaceExplicitAlgorithm <: ODEInterfaceAlgorithm end
"""
    dopri5()

Hairer's classic implementation of the Dormand-Prince 4/5 method.

An explicit Runge-Kutta method of order 5(4) for non-stiff problems. This is one of the most
widely used explicit methods for solving ODEs and is generally efficient for problems where
evaluations of the right-hand side are relatively cheap.

## Solver Properties
- Order: 5
- Adaptive: Yes
- Type: Explicit Runge-Kutta
- Suitable for: Non-stiff problems
"""
struct dopri5 <: ODEInterfaceExplicitAlgorithm end
SciMLBase.alg_order(alg::dopri5) = 5

"""
    dop853()

Explicit Runge-Kutta 8(5,3) method by Dormand-Prince.

A high-order explicit Runge-Kutta method for non-stiff problems requiring high accuracy.
Uses adaptive step size control and is particularly effective when tight tolerances are needed.

## Solver Properties
- Order: 8
- Adaptive: Yes
- Type: Explicit Runge-Kutta
- Suitable for: Non-stiff problems requiring high accuracy
"""
struct dop853 <: ODEInterfaceExplicitAlgorithm end
SciMLBase.alg_order(alg::dop853) = 8

"""
    odex()

GBS (Gragg-Bulirsch-Stoer) extrapolation algorithm based on the midpoint rule.

An extrapolation method that can achieve very high accuracy for smooth non-stiff problems.
Particularly efficient for problems where high accuracy is required and the solution is smooth.

## Solver Properties
- Order: up to 12
- Adaptive: Yes
- Type: Extrapolation
- Suitable for: Smooth non-stiff problems requiring high accuracy
"""
struct odex <: ODEInterfaceExplicitAlgorithm end
SciMLBase.alg_order(alg::odex) = 12

"""
    seulex(; jac_lower = nothing, jac_upper = nothing)

Extrapolation algorithm based on the linearly implicit Euler method.

An implicit extrapolation method suitable for stiff and differential-algebraic problems.
Can handle mildly stiff to stiff systems efficiently with automatic stiffness detection.

## Arguments
- `jac_lower`: Lower bandwidth of the Jacobian for banded matrices
- `jac_upper`: Upper bandwidth of the Jacobian for banded matrices

## Solver Properties
- Order: up to 12
- Adaptive: Yes
- Type: Implicit extrapolation
- Suitable for: Stiff problems and DAEs
"""
struct seulex{T} <: ODEInterfaceImplicitAlgorithm
    jac_lower::T
    jac_upper::T
end
SciMLBase.alg_order(alg::seulex) = 12

"""
    radau(; jac_lower = nothing, jac_upper = nothing)

Implicit Runge-Kutta (Radau IIA) method of variable order between 5 and 13.

A high-order implicit Runge-Kutta method particularly well-suited for stiff ODEs and DAEs.
The Radau IIA methods have excellent stability properties and can handle very stiff problems
efficiently.

## Arguments
- `jac_lower`: Lower bandwidth of the Jacobian for banded matrices
- `jac_upper`: Upper bandwidth of the Jacobian for banded matrices

## Solver Properties
- Order: 5 to 13 (variable)
- Adaptive: Yes
- Type: Implicit Runge-Kutta (Radau IIA)
- Suitable for: Stiff problems and DAEs
"""
struct radau{T} <: ODEInterfaceImplicitAlgorithm
    jac_lower::T
    jac_upper::T
end
SciMLBase.alg_order(alg::radau) = 13

"""
    radau5(; jac_lower = nothing, jac_upper = nothing)

Implicit Runge-Kutta method (Radau IIA) of order 5.

A robust implicit method for stiff ODEs and DAEs. This is the fixed order 5 version of the
Radau method, often preferred when a specific order is desired or for comparison purposes.

## Arguments
- `jac_lower`: Lower bandwidth of the Jacobian for banded matrices
- `jac_upper`: Upper bandwidth of the Jacobian for banded matrices

## Solver Properties
- Order: 5
- Adaptive: Yes
- Type: Implicit Runge-Kutta (Radau IIA)
- Suitable for: Stiff problems and DAEs
"""
struct radau5{T} <: ODEInterfaceImplicitAlgorithm
    jac_lower::T
    jac_upper::T
end
SciMLBase.alg_order(alg::radau5) = 5

"""
    rodas(; jac_lower = nothing, jac_upper = nothing)

Rosenbrock method of order 4(3).

A Rosenbrock-type method for stiff ODEs and DAEs. Rosenbrock methods are semi-implicit,
using only one LU decomposition per step, which can make them more efficient than fully
implicit methods for moderately stiff problems.

## Arguments
- `jac_lower`: Lower bandwidth of the Jacobian for banded matrices
- `jac_upper`: Upper bandwidth of the Jacobian for banded matrices

## Solver Properties
- Order: 4
- Adaptive: Yes
- Type: Rosenbrock
- Suitable for: Stiff problems and DAEs
"""
struct rodas{T} <: ODEInterfaceImplicitAlgorithm
    jac_lower::T
    jac_upper::T
end
SciMLBase.alg_order(alg::rodas) = 4

"""
    ddeabm()

Adams-Bashforth-Moulton predictor-corrector method with variable order (1 to 12).

A multistep method based on Adams formulas, using a predictor-corrector approach.
Efficient for smooth non-stiff problems, particularly when many function evaluations
are expensive, as it reuses previous function evaluations.

## Solver Properties
- Order: 1 to 12 (variable)
- Adaptive: Yes
- Type: Adams-Bashforth-Moulton multistep
- Suitable for: Smooth non-stiff problems
"""
struct ddeabm <: ODEInterfaceExplicitAlgorithm end
SciMLBase.alg_order(alg::ddeabm) = 12

"""
    ddebdf(; jac_lower = nothing, jac_upper = nothing)

Backward Differentiation Formula (BDF) with variable order (1 to 5).

A multistep method particularly effective for stiff ODEs and DAEs. BDF methods are
the industry standard for stiff problems in many scientific computing applications.

## Arguments
- `jac_lower`: Lower bandwidth of the Jacobian for banded matrices
- `jac_upper`: Upper bandwidth of the Jacobian for banded matrices

## Solver Properties
- Order: 1 to 5 (variable)
- Adaptive: Yes
- Type: BDF multistep
- Suitable for: Stiff problems and DAEs
"""
struct ddebdf{T} <: ODEInterfaceImplicitAlgorithm
    jac_lower::T
    jac_upper::T
end
SciMLBase.alg_order(alg::ddebdf) = 5

seulex(; jac_lower = nothing, jac_upper = nothing) = seulex(jac_lower, jac_upper)
radau(; jac_lower = nothing, jac_upper = nothing) = radau(jac_lower, jac_upper)
radau5(; jac_lower = nothing, jac_upper = nothing) = radau5(jac_lower, jac_upper)
rodas(; jac_lower = nothing, jac_upper = nothing) = rodas(jac_lower, jac_upper)
ddebdf(; jac_lower = nothing, jac_upper = nothing) = ddebdf(jac_lower, jac_upper)
