__precompile__()

module ODEInterfaceDiffEq
  using DiffEqBase, ODEInterface, Compat

  import DiffEqBase: solve

  @compat const KW = Dict{Symbol,Any}

  include("algorithms.jl")
  include("solve.jl")

  export ODEInterfaceAlgorithm, dopri5, dop853, odex, seulex, radau, radau5, rodas


end # module
