using ODEInterfaceDiffEq, BenchmarkTools

const SUITE = BenchmarkGroup()

function lorenz!(du, u, p, t)
    du[1] = 10.0 * (u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
    return nothing
end
prob = ODEProblem(lorenz!, [1.0, 0.0, 0.0], (0.0, 50.0))

# Stiff problem
function rober!(du, u, p, t)
    du[1] = -0.04u[1] + 1.0e4 * u[2] * u[3]
    du[2] = 0.04u[1] - 1.0e4 * u[2] * u[3] - 3.0e7 * u[2]^2
    du[3] = 3.0e7 * u[2]^2
    return nothing
end
prob_stiff = ODEProblem(rober!, [1.0, 0.0, 0.0], (0.0, 1.0e5))

# =============================================================================
# Non-stiff solvers
# =============================================================================

SUITE["nonstiff"] = BenchmarkGroup()

SUITE["nonstiff"]["dopri5"] = @benchmarkable solve($prob, dopri5())
SUITE["nonstiff"]["dop853"] = @benchmarkable solve($prob, dop853())
SUITE["nonstiff"]["odex"] = @benchmarkable solve($prob, odex())

# =============================================================================
# Stiff solvers
# =============================================================================

SUITE["stiff"] = BenchmarkGroup()

SUITE["stiff"]["radau5"] = @benchmarkable solve(
    $prob_stiff, radau5(); saveat = 1.0e3
)
SUITE["stiff"]["rodas"] = @benchmarkable solve(
    $prob_stiff, rodas(); saveat = 1.0e3
)
SUITE["stiff"]["seulex"] = @benchmarkable solve(
    $prob_stiff, seulex(); saveat = 1.0e3
)
