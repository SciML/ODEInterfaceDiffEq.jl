# DAE Initialization support for ODEInterfaceDiffEq
# Following the pattern from Sundials.jl:
# https://github.com/SciML/Sundials.jl/blob/master/src/common_interface/initialize.jl

import SciMLBase: OverrideInit, NoInit, CheckInit, has_initialization_data
import DiffEqBase: DefaultInit

# Re-export initialization algorithms (including DefaultInit from DiffEqBase)
export OverrideInit, NoInit, CheckInit, DefaultInit

# DefaultInit: OverrideInit → CheckInit pattern (matching Sundials v5)
# First run OverrideInit to compute consistent initial conditions,
# then run CheckInit to verify the algebraic constraints are satisfied.
function DiffEqBase.initialize_dae!(
        integrator::ODEInterfaceIntegrator,
        initializealg::DefaultInit
    )
    prob = integrator.sol.prob

    # First: OverrideInit to compute consistent initial conditions
    if has_initialization_data(prob.f)
        DiffEqBase.initialize_dae!(integrator, OverrideInit())

        # Check if OverrideInit failed
        if integrator.sol.retcode == ReturnCode.InitialFailure
            return nothing
        end
    end

    # Then: CheckInit to verify algebraic constraints are satisfied
    DiffEqBase.initialize_dae!(integrator, CheckInit())

    return nothing
end

# NoInit: Do nothing, assume initial conditions are correct
function DiffEqBase.initialize_dae!(
        integrator::ODEInterfaceIntegrator,
        initializealg::NoInit
    )
    # No-op: initial conditions are assumed to be correct
    return nothing
end

# CheckInit: Verify that initial conditions satisfy the algebraic constraints
function DiffEqBase.initialize_dae!(
        integrator::ODEInterfaceIntegrator,
        initializealg::CheckInit
    )
    prob = integrator.sol.prob
    f = prob.f
    M = f.mass_matrix

    # If no mass matrix or identity, no algebraic constraints to check
    M == I && return nothing

    u0 = integrator.u
    p = integrator.p
    t = integrator.t

    # Find algebraic equations (rows of M that are all zeros)
    algebraic_eqs = [all(iszero, M[i, :]) for i in axes(M, 1)]

    # If no algebraic equations, nothing to check
    !any(algebraic_eqs) && return nothing

    # Evaluate the RHS
    tmp = similar(u0)
    f(tmp, u0, p, t)

    # Check residuals of algebraic equations only
    abstol = integrator.sol.prob isa DiffEqBase.AbstractODEProblem ?
             get(Dict(integrator.opts.callback.discrete_callbacks), :abstol, 1e-8) : 1e-8

    # Try to get abstol from the solve options, fallback to default
    # Note: ODEInterface doesn't expose tolerances through the integrator in the same way
    abstol = 1e-8  # Default tolerance

    max_residual = maximum(abs.(tmp[algebraic_eqs]))

    if max_residual > abstol
        error(
            """
            DAE initialization failed with CheckInit: Initial conditions do not satisfy the algebraic constraints.

            The maximum residual in algebraic equations is $(max_residual), which exceeds the tolerance $(abstol).

            To resolve this issue, you have several options:
            1. Fix your initial conditions to satisfy the algebraic constraints (M * du/dt = f(u, p, t), where algebraic rows have M[i,:] = 0)
            2. If using ModelingToolkit, use: initializealg = OverrideInit()
            3. Use initializealg = NoInit() to skip initialization checks (use with caution)

            Example:
            solve(prob, radau5(); initializealg = OverrideInit())
            """
        )
    end

    return nothing
end

# OverrideInit: Use SciMLBase's initialization system (e.g., from ModelingToolkit)
function DiffEqBase.initialize_dae!(
        integrator::ODEInterfaceIntegrator,
        initializealg::OverrideInit
    )
    prob = integrator.sol.prob
    f = prob.f

    # If no initialization data, nothing to do
    if !has_initialization_data(f)
        return nothing
    end

    # Get initial values using SciMLBase's initialization system
    # Note: ODEInterface integrators don't have a built-in nonlinear solver,
    # so we rely on the user providing one via nlsolve_alg in OverrideInit,
    # or the initialization being trivial (no solve needed)
    isinplace = Val(SciMLBase.isinplace(prob))

    # Default tolerances
    abstol = 1e-8
    reltol = 1e-8

    u0, p, success = SciMLBase.get_initial_values(
        prob, integrator, f, initializealg, isinplace;
        abstol = abstol, reltol = reltol
    )

    if !success
        integrator.sol = SciMLBase.solution_new_retcode(
            integrator.sol,
            ReturnCode.InitialFailure
        )
        return nothing
    end

    # Update integrator state
    if SciMLBase.isinplace(prob)
        integrator.u .= u0
        if length(integrator.sol.u) >= 1 && !isempty(integrator.sol.u)
            integrator.sol.u[1] .= u0
        end
    else
        # For out-of-place problems, we need to handle this differently
        # since integrator.u might be immutable
        integrator.u .= u0
        if length(integrator.sol.u) >= 1 && !isempty(integrator.sol.u)
            integrator.sol.u[1] = u0
        end
    end

    # Update parameters if they changed
    if p !== integrator.p
        # Note: Updating parameters in the integrator for ODEInterface
        # is tricky because the parameters are captured in the closure.
        # This is a known limitation - the parameters in the function closure
        # won't be updated. For full support, the problem would need to be re-created.
        @warn "OverrideInit updated parameters, but ODEInterface solvers may not reflect parameter changes during integration. Consider re-creating the problem with the new parameters."
    end

    integrator.u_modified = true

    return nothing
end
