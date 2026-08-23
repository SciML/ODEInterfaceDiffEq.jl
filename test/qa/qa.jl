using ODEInterfaceDiffEq, SciMLTesting, Test

# The SciML common interface ODEInterfaceDiffEq deliberately reexports so that
# `using ODEInterfaceDiffEq` is enough to build and solve an ODE problem, as the README
# example does. Owned and documented upstream; kept in sync with the reexport `export`
# block in src/ODEInterfaceDiffEq.jl.
const REEXPORTS = (
    :CallbackSet, :CheckInit, :ContinuousCallback, :DECallback, :DEIntegrator,
    :DEStats, :DefaultInit, :DiscreteCallback, :EnsembleAnalysis, :EnsembleDistributed,
    :EnsembleProblem, :EnsembleSerial, :EnsembleSolution, :EnsembleSplitThreads,
    :EnsembleSummary, :EnsembleThreads, :NoInit, :NullParameters, :ODEFunction,
    :ODEProblem, :ODESolution, :OverrideInit, :ReturnCode, :VectorContinuousCallback,
    :change_t_via_interpolation!, :get_tmp_cache, :init, :remake, :savevalues!,
    :set_proposed_dt!, :solve, :solve!, :step!, :successful_retcode, :terminate!,
    :u_modified!,
)

run_qa(ODEInterfaceDiffEq; reexports_allow = REEXPORTS)

@testset "Reexport surface" begin
    # Every approved reexport must actually be reachable from `using ODEInterfaceDiffEq`,
    # so the allow-list cannot drift into approving names the package no longer provides.
    @testset "$name" for name in REEXPORTS
        @test name in names(ODEInterfaceDiffEq)
        @test isdefined(@__MODULE__, name)
    end
end
