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
    using SciMLBase: CallbackSet, ReturnCode, VectorContinuousCallback, check_keywords,
        warn_compat
    using SciMLLogging: SciMLLogging, @SciMLMessage
    import DiffEqBase: initialize!, savevalues!

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

end # module
