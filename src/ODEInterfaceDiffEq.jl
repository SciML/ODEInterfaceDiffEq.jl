__precompile__()

module ODEInterfaceDiffEq

    import DiffEqBase

    import Compat
    import FunctionWrappers
    import ODEInterface
    import SciMLBase
    import SciMLStructures
    using DataStructures: BinaryMaxHeap, BinaryMinHeap, counter
    using LinearAlgebra: I
    using SciMLBase: check_keywords, warn_compat
    using SciMLLogging: SciMLLogging, @SciMLMessage
    import DiffEqBase: initialize!, savevalues!

    # The SciML common interface that ODEInterfaceDiffEq reexports (see the second
    # `export` below), so that `using ODEInterfaceDiffEq` on its own is enough to build
    # an ODE problem, solve it with one of the ODEInterface solvers, drive the
    # integrator from a callback, and inspect the result -- which is exactly what the
    # README example does. The integrator functions listed here are the ones this
    # package implements methods for. Every name stays owned and documented upstream.
    using SciMLBase: CallbackSet, CheckInit, ContinuousCallback, DECallback, DEIntegrator,
        DEStats, DiscreteCallback, EnsembleAnalysis, EnsembleDistributed, EnsembleProblem,
        EnsembleSerial, EnsembleSolution, EnsembleSplitThreads, EnsembleSummary,
        EnsembleThreads, NoInit, NullParameters, ODEFunction, ODEProblem, ODESolution,
        OverrideInit, ReturnCode, VectorContinuousCallback, change_t_via_interpolation!,
        get_tmp_cache, init, remake, set_proposed_dt!, solve, solve!, step!,
        successful_retcode, terminate!, u_modified!
    using DiffEqBase: DefaultInit

    const warnkeywords = (
        :save_idxs, :d_discontinuities, :unstable_check, :tstops,
        :calck, :progress, :dense, :save_start,
    )

    function __init__()
        return global warnlist = Set(warnkeywords)
    end

    const KW = Dict{Symbol, Any}

    include("algorithms.jl")
    include("integrator_types.jl")
    include("integrator_utils.jl")
    include("solve.jl")
    include("initialize.jl")

    using PrecompileTools: @compile_workload, @setup_workload

    @setup_workload begin
        @compile_workload begin
            dopri5()
            dop853()
            radau5()
            rodas()
            ddeabm()
            ddebdf()
        end
    end

    export ODEInterfaceAlgorithm, dopri5, dop853, odex, seulex, radau, radau5, rodas,
        ddeabm, ddebdf

    # Reexported SciML common interface; approved via `reexports_allow` in test/qa/qa.jl.
    export CallbackSet, CheckInit, ContinuousCallback, DECallback, DEIntegrator,
        DEStats, DefaultInit, DiscreteCallback, EnsembleAnalysis, EnsembleDistributed,
        EnsembleProblem, EnsembleSerial, EnsembleSolution, EnsembleSplitThreads,
        EnsembleSummary, EnsembleThreads, NoInit, NullParameters, ODEFunction, ODEProblem,
        ODESolution, OverrideInit, ReturnCode, VectorContinuousCallback,
        change_t_via_interpolation!, get_tmp_cache, init, remake, savevalues!,
        set_proposed_dt!, solve, solve!, step!, successful_retcode, terminate!,
        u_modified!

end # module
