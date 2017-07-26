__precompile__()

module ODEInterfaceDiffEq

using DiffEqBase, ODEInterface, Compat

import DiffEqBase: solve

const warnkeywords =
    (:save_idxs, :d_discontinuities, :unstable_check,
     :calck, :progress, :timeseries_steps, :dense)

function __init__()
    const global warnlist = Set(warnkeywords)
end

@compat const KW = Dict{Symbol,Any}

include("algorithms.jl")
include("solve.jl")

export ODEInterfaceAlgorithm, dopri5, dop853, odex, seulex, radau, radau5, rodas,
       ddeabm, ddebdf

end # module
