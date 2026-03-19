using ODEInterfaceDiffEq, DiffEqBase, SciMLBase
using Test
using LinearAlgebra

@testset "DAE Initialization" begin
    # Simple Robertson chemical kinetics model as a DAE (with mass matrix)
    # dy1/dt = -0.04*y1 + 1e4*y2*y3
    # dy2/dt = 0.04*y1 - 1e4*y2*y3 - 3e7*y2^2
    # 0 = y1 + y2 + y3 - 1 (conservation constraint)

    function robertson!(du, u, p, t)
        y1, y2, y3 = u
        du[1] = -0.04 * y1 + 1e4 * y2 * y3
        du[2] = 0.04 * y1 - 1e4 * y2 * y3 - 3e7 * y2^2
        du[3] = y1 + y2 + y3 - 1.0  # algebraic equation
        return nothing
    end

    # Mass matrix: [1 0 0; 0 1 0; 0 0 0]
    M = [1.0 0 0; 0 1 0; 0 0 0]

    @testset "NoInit - skip initialization" begin
        # With correct initial conditions
        u0_correct = [1.0, 0.0, 0.0]  # satisfies y1 + y2 + y3 = 1
        f = ODEFunction(robertson!; mass_matrix = M)
        prob = ODEProblem(f, u0_correct, (0.0, 1e-3))

        sol = solve(prob, radau5(); initializealg = NoInit())
        @test sol.retcode == ReturnCode.Success

        # With incorrect initial conditions - NoInit should still work (no check)
        u0_wrong = [1.0, 1.0, 1.0]  # does NOT satisfy y1 + y2 + y3 = 1
        prob_wrong = ODEProblem(f, u0_wrong, (0.0, 1e-3))

        # NoInit should not throw, even with wrong ICs
        sol_wrong = solve(prob_wrong, radau5(); initializealg = NoInit())
        # The solver might fail, but not during initialization
        @test true  # Just checking NoInit doesn't throw
    end

    @testset "CheckInit - verify constraints" begin
        u0_correct = [1.0, 0.0, 0.0]  # satisfies y1 + y2 + y3 = 1
        f = ODEFunction(robertson!; mass_matrix = M)
        prob = ODEProblem(f, u0_correct, (0.0, 1e-3))

        sol = solve(prob, radau5(); initializealg = CheckInit())
        @test sol.retcode == ReturnCode.Success

        # With incorrect initial conditions
        u0_wrong = [1.0, 1.0, 1.0]  # does NOT satisfy y1 + y2 + y3 = 1
        prob_wrong = ODEProblem(f, u0_wrong, (0.0, 1e-3))

        @test_throws ErrorException solve(prob_wrong, radau5(); initializealg = CheckInit())
    end

    @testset "DefaultInit - dispatch based on problem" begin
        # Without initialization_data, should use CheckInit
        u0_correct = [1.0, 0.0, 0.0]
        f = ODEFunction(robertson!; mass_matrix = M)
        prob = ODEProblem(f, u0_correct, (0.0, 1e-3))

        sol = solve(prob, radau5(); initializealg = DefaultInit())
        @test sol.retcode == ReturnCode.Success
    end

    @testset "No mass matrix - identity case" begin
        # Simple ODE without mass matrix (not a DAE)
        function simple_ode!(du, u, p, t)
            du[1] = -u[1]
            return nothing
        end

        u0 = [1.0]
        f = ODEFunction(simple_ode!)
        prob = ODEProblem(f, u0, (0.0, 1.0))

        # All initialization algorithms should work
        sol_noinit = solve(prob, radau5(); initializealg = NoInit())
        @test sol_noinit.retcode == ReturnCode.Success

        sol_check = solve(prob, radau5(); initializealg = CheckInit())
        @test sol_check.retcode == ReturnCode.Success

        sol_default = solve(prob, radau5(); initializealg = DefaultInit())
        @test sol_default.retcode == ReturnCode.Success
    end

    @testset "Multiple solvers with initialization" begin
        u0_correct = [1.0, 0.0, 0.0]
        f = ODEFunction(robertson!; mass_matrix = M)
        prob = ODEProblem(f, u0_correct, (0.0, 1e-3))

        # Test with different implicit solvers that support mass matrices
        for alg in [radau5(), radau(), rodas(), seulex(), ddebdf()]
            sol = solve(prob, alg; initializealg = CheckInit())
            @test sol.retcode == ReturnCode.Success
        end
    end
end
