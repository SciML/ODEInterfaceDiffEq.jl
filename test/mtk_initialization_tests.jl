# ModelingToolkit initialization tests for ODEInterfaceDiffEq
# Based on Sundials.jl/test/common_interface/initialization.jl

using ODEInterfaceDiffEq, DiffEqBase, SciMLBase, Test
using ModelingToolkit
using NonlinearSolve
using SymbolicIndexingInterface
using ModelingToolkit: t_nounits as t, D_nounits as D

@testset "MTK ODE Initialization" begin
    # ODE with missing parameters that need to be determined via initialization
    # System: dx/dt = p*y + q*t, dy/dt = 5x + q
    # With initialization equations to find p and q from initial state constraints
    @variables x(t) [guess = 1.0] y(t) [guess = 1.0]
    @parameters p = missing [guess = 1.0] q = missing [guess = 1.0]

    @mtkcompile sys = System(
        [D(x) ~ p * y + q * t, D(y) ~ 5x + q],
        t;
        initialization_eqs = [p^2 + q^2 ~ 3, x^3 + y^3 ~ 5]
    )

    @testset "IIP: $iip" for iip in [true, false]
        prob = ODEProblem{iip}(sys, [x => 1.0, p => 1.0], (0.0, 0.1))

        # Test with implicit solvers that support OverrideInit
        @testset "$alg" for alg in [radau5, radau, rodas, seulex]
            sol = solve(prob, alg())
            @test SciMLBase.successful_retcode(sol)
            @test sol[x, 1] ≈ 1.0
            @test sol[y, 1] ≈ cbrt(4)
            @test sol.ps[p] ≈ 1.0
            @test sol.ps[q] ≈ sqrt(2)
        end
    end
end

@testset "MTK DAE Initialization (Mass Matrix)" begin
    # DAE system using mass matrix formulation
    # System: dx/dt = p*y + q*t, x^3 + y^3 = 5 (algebraic constraint)
    # With initialization equations to find the missing parameter q
    @variables x(t) [guess = 1.0] y(t) [guess = 1.0]
    @parameters p = missing [guess = 1.0] q = missing [guess = 1.0]

    @mtkcompile sys = System(
        [D(x) ~ p * y + q * t, x^3 + y^3 ~ 5],
        t;
        initialization_eqs = [p^2 + q^2 ~ 3]
    )

    @testset "OverrideInit" begin
        # Provide D(x) guess and p value - initialization should determine x, y, q
        prob = ODEProblem(sys, [D(x) => cbrt(4), p => 1.0], (0.0, 0.1))

        @testset "$alg" for alg in [radau5, radau]
            sol = solve(prob, alg())
            @test SciMLBase.successful_retcode(sol)
            @test sol[x, 1] ≈ 1.0
            @test sol[y, 1] ≈ cbrt(4)
            @test sol.ps[p] ≈ 1.0
            @test sol.ps[q] ≈ sqrt(2)
        end
    end

    @testset "CheckInit" begin
        prob = ODEProblem(sys, [D(x) => cbrt(4), p => 1.0], (0.0, 0.1))

        # CheckInit should fail because the MTK problem has initialization_data
        # and the user-provided values need processing via OverrideInit
        @test_throws Any solve(prob, radau5(); initializealg = SciMLBase.CheckInit())
    end
end

@testset "MTK Pendulum DAE" begin
    # Classic pendulum as index-3 DAE (reduced to index-1 via ModelingToolkit)
    # This tests a more complex DAE initialization scenario
    @parameters g = 9.81 L = 1.0
    @variables begin
        x(t), [guess = 1.0]
        y(t), [guess = 0.0]
        vx(t), [guess = 0.0]
        vy(t), [guess = 0.0]
        λ(t), [guess = 0.0]  # Lagrange multiplier
    end

    eqs = [
        D(x) ~ vx
        D(y) ~ vy
        D(vx) ~ -2λ * x
        D(vy) ~ -2λ * y - g
        x^2 + y^2 ~ L^2  # algebraic constraint
    ]

    @mtkcompile pend = System(eqs, t)

    # Initial condition: pendulum at angle θ₀ = π/6 from vertical
    θ₀ = π / 6
    x0 = L * sin(θ₀)
    y0 = -L * cos(θ₀)

    prob = ODEProblem(pend, [x => x0, y => y0, vx => 0.0, vy => 0.0], (0.0, 0.5))

    @testset "$alg" for alg in [radau5, radau]
        sol = solve(prob, alg())
        @test SciMLBase.successful_retcode(sol)
        # Check that the constraint is satisfied throughout
        @test all(abs.((sol[x] .^ 2 .+ sol[y] .^ 2) .- L^2) .< 1.0e-4)
    end
end
