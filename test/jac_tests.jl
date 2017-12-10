using ODEInterfaceDiffEq, DiffEqProblemLibrary, DiffEqBase
using Base.Test

jac_called = false

function Lotka(t,u,du)
  du[1] = u[1] - u[1] * u[2] # REPL[7], line 3:
  du[2] = -3 * u[2] + 1 * u[1] * u[2]
  nothing
end

function Lotka(::Type{Val{:jac}},t,u,J)
  global jac_called
  jac_called = true
  J[1,1] = 1.0 - u[2]
  J[1,2] = -u[1]
  J[2,1] = 1 * u[2]
  J[2,2] = -3 + u[1]
  nothing
end

prob = ODEProblem(Lotka,ones(2),(0.0,2.0))

sol =solve(prob,radau5();dt=1//2^(4))

@test jac_called == true
