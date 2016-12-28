module ODEInterfaceDiffEq
  using DiffEqBase, ODEInterface

  import DiffEqBase: solve

  typealias KW Dict{Symbol,Any}

  include("algorithms.jl")
  include("solve.jl")

  export ODEInterfaceAlgorithm, dopri5, dop853, odex, seulex, radau, radau5, rodas


end # module
