# ODEInterfaceDiffEq

[![Join the chat at https://gitter.im/JuliaDiffEq/Lobby](https://badges.gitter.im/JuliaDiffEq/Lobby.svg)](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Build Status](https://travis-ci.org/JuliaDiffEq/ODEInterfaceDiffEq.jl.svg?branch=master)](https://travis-ci.org/JuliaDiffEq/ODEInterfaceDiffEq.jl)
[![Coverage Status](https://coveralls.io/repos/JuliaDiffEq/ODEInterfaceDiffEq.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/JuliaDiffEq/ODEInterfaceDiffEq.jl?branch=master)
[![codecov.io](http://codecov.io/github/JuliaDiffEq/ODEInterfaceDiffEq.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaDiffEq/ODEInterfaceDiffEq.jl?branch=master)

[![ODEInterfaceDiffEq](http://pkg.julialang.org/badges/ODEInterfaceDiffEq_0.5.svg)](http://pkg.julialang.org/?pkg=ODEInterfaceDiffEq)
[![ODEInterfaceDiffEq](http://pkg.julialang.org/badges/ODEInterfaceDiffEq_0.6.svg)](http://pkg.julialang.org/?pkg=ODEInterfaceDiffEq)

This package contains bindings for ODEInterface.jl to allow it to be used with the
JuliaDiffEq common interface. For more information on using the solvers from this
package, see the [DifferentialEquations.jl documentation](https://juliadiffeq.github.io/DiffEqDocs.jl/latest/).

## Installation

A standard installation on MacOSX and Linux should work. On Windows, you need to install mingw32 compilers and add them to the path. [MingW32 can be found here](http://www.mingw.org/). Then add the path to your environment variables. An example path is:

```
C:\Program Files\mingw-w64\x86_64-6.1.0-posix-seh-rt_v5-rev0\mingw64\bin
```
