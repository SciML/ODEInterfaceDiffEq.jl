# ODEInterfaceDiffEq.jl

```@docs
ODEInterfaceAlgorithm
dopri5
dop853
odex
seulex
radau
radau5
rodas
ddeabm
ddebdf
```

## Reexported SciML common interface

`using ODEInterfaceDiffEq` also brings in the parts of the SciML common interface needed
to build a problem, solve it with one of the solvers above, and inspect the result, so
they do not have to be imported separately. These names are owned and documented by
[SciMLBase](https://docs.sciml.ai/SciMLBase/stable/) and
[DiffEqBase](https://docs.sciml.ai/DiffEqDocs/stable/) — `DefaultInit` comes from
DiffEqBase, everything else from SciMLBase:

  - Problems: `ODEProblem`, `EnsembleProblem`
  - Functions: `ODEFunction`
  - Solutions: `ODESolution`, `EnsembleSolution`, `EnsembleSummary`, `DEStats`
  - Ensemble algorithms: `EnsembleSerial`, `EnsembleThreads`, `EnsembleDistributed`,
    `EnsembleSplitThreads`, and the `EnsembleAnalysis` module
  - Solving: `solve`, `solve!`, `init`, `step!`, `remake`
  - Integrator interface: `DEIntegrator`, `savevalues!`, `change_t_via_interpolation!`,
    `get_tmp_cache`, `set_proposed_dt!`, `u_modified!`, `terminate!`
  - Return status: `ReturnCode`, `successful_retcode`
  - Callbacks: `ContinuousCallback`, `DiscreteCallback`, `VectorContinuousCallback`,
    `CallbackSet`, `DECallback`
  - Initialization algorithms: `DefaultInit`, `NoInit`, `CheckInit`, `OverrideInit`
  - `NullParameters`

The integrator functions listed above are exactly the ones ODEInterfaceDiffEq implements
methods for, and the initialization algorithms are exactly the ones it implements
`initialize_dae!` for. `BrownFullBasicInit` and `ShampineCollocationInit` are
deliberately not reexported, since they have no method here. Anything else from SciMLBase
or DiffEqBase must be imported from that package directly.
