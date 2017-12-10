__precompile__()

module ODEInterfaceDiffEq

using Reexport
@reexport using DiffEqBase

using ODEInterface, Compat, DataStructures

import DiffEqBase: solve

const warnkeywords =
    (:save_idxs, :d_discontinuities, :unstable_check, :tstops,
     :calck, :progress, :timeseries_steps, :dense,:save_start)

function __init__()
    const global warnlist = Set(warnkeywords)
end

@compat const KW = Dict{Symbol,Any}

include("algorithms.jl")
include("integrator_types.jl")
include("integrator_utils.jl")
include("callbacks.jl")
include("solve.jl")

export ODEInterfaceAlgorithm, dopri5, dop853, odex, seulex, radau, radau5, rodas,
       ddeabm, ddebdf

end # module
