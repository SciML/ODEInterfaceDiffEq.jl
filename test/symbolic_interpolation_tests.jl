using ODEInterfaceDiffEq, Test
using SciMLBase: SymbolCache

# The ODEInterface integrator only exists inside a callback, so the symbolic indexing is
# probed from an `affect!`. `SymbolCache` resolves an `Expr` as an observed quantity, so
# `:(x + y)` is not a component of the state vector.
f = ODEFunction(
    (du, u, p, t) -> (du[1] = u[2]; du[2] = -u[1]);
    sys = SymbolCache([:x, :y], Symbol[], :t)
)
prob = ODEProblem(f, [1.0, 0.0], (0.0, 10.0))

@testset "Symbolic idxs inside a callback" begin
    probed = Ref(false)
    affect! = function (integrator)
        t = integrator.t
        # A symbolic index evaluates the whole state and then projects, while an integer
        # index reads that component out of the raw output, so these agree to rounding.
        @test integrator(t; idxs = :x) ≈ integrator(t; idxs = 1)
        @test integrator(t; idxs = :y) ≈ integrator(t; idxs = 2)
        @test integrator(t; idxs = :(x + y)) ≈ sum(integrator(t))
        @test integrator(t; idxs = [:x, :y]) ≈ integrator(t)

        # Symbolic routing does not disturb the existing integer and `nothing` paths
        @test integrator(t) == integrator(t; idxs = nothing)
        @test integrator(t; idxs = [1, 2]) ≈ integrator(t)
        probed[] = true
        return nothing
    end

    # A ContinuousCallback is what puts the solver in dense-output mode, which the
    # interpolant this test exercises requires.
    solve(
        prob, dopri5();
        callback = ContinuousCallback((u, t, integrator) -> u[1], affect!),
        dtmax = 0.5
    )
    @test probed[]
end
