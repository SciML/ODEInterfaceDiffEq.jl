# ODEInterface.jl Algorithms

@compat abstract type ODEInterfaceAlgorithm <: AbstractODEAlgorithm end
immutable dopri5 <: ODEInterfaceAlgorithm end
immutable dop853 <: ODEInterfaceAlgorithm end
immutable odex <: ODEInterfaceAlgorithm end
immutable seulex <: ODEInterfaceAlgorithm end
immutable radau <: ODEInterfaceAlgorithm end
immutable radau5 <: ODEInterfaceAlgorithm end
immutable rodas <: ODEInterfaceAlgorithm end
