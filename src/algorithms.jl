# ODEInterface.jl Algorithms

abstract type ODEInterfaceAlgorithm <: DiffEqBase.AbstractODEAlgorithm end
struct dopri5 <: ODEInterfaceAlgorithm end
struct dop853 <: ODEInterfaceAlgorithm end
struct odex <: ODEInterfaceAlgorithm end
struct seulex <: ODEInterfaceAlgorithm end
struct radau <: ODEInterfaceAlgorithm end
struct radau5 <: ODEInterfaceAlgorithm end
struct rodas <: ODEInterfaceAlgorithm end
struct ddeabm <: ODEInterfaceAlgorithm end
struct ddebdf <: ODEInterfaceAlgorithm end
