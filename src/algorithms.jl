# ODEInterface.jl Algorithms

abstract type ODEInterfaceAlgorithm <: DiffEqBase.AbstractODEAlgorithm end
abstract type ODEInterfaceImplicitAlgorithm <: ODEInterfaceAlgorithm end
abstract type ODEInterfaceExplicitAlgorithm <: ODEInterfaceAlgorithm end
"""
dopri5: Hairer's classic implementation of the Dormand-Prince 4/5 method.
"""
struct dopri5 <: ODEInterfaceExplicitAlgorithm end
SciMLBase.alg_order(alg::dopri5) = 5

"""
dop853: Explicit Runge-Kutta 8(5,3) by Dormand-Prince.
"""
struct dop853 <: ODEInterfaceExplicitAlgorithm end
SciMLBase.alg_order(alg::dop853) = 8

"""
odex: GBS extrapolation-algorithm based on the midpoint rule.
"""
struct odex <: ODEInterfaceExplicitAlgorithm end
SciMLBase.alg_order(alg::odex) = 12

"""
seulex: Extrapolation-algorithm based on the linear implicit Euler method.
"""
struct seulex{T} <: ODEInterfaceImplicitAlgorithm
    jac_lower::T
    jac_upper::T
end
SciMLBase.alg_order(alg::seulex) = 12

"""
radau: Implicit Runge-Kutta (Radau IIA) of variable order between 5 and 13.
"""
struct radau{T} <: ODEInterfaceImplicitAlgorithm
    jac_lower::T
    jac_upper::T
end
SciMLBase.alg_order(alg::radau) = 13

"""
radau5: Implicit Runge-Kutta method (Radau IIA) of order 5.
"""
struct radau5{T} <: ODEInterfaceImplicitAlgorithm
    jac_lower::T
    jac_upper::T
end
SciMLBase.alg_order(alg::radau5) = 5

"""
rodas: Rosenbrock 4(3) method.
"""
struct rodas{T} <: ODEInterfaceImplicitAlgorithm
    jac_lower::T
    jac_upper::T
end
SciMLBase.alg_order(alg::rodas) = 4

"""
ddeabm: Adams-Bashforth-Moulton Predictor-Corrector method (order between 1 and 12)
"""
struct ddeabm <: ODEInterfaceExplicitAlgorithm end
SciMLBase.alg_order(alg::ddeabm) = 12

"""
ddebdf: Backward Differentiation Formula (orders between 1 and 5)
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
