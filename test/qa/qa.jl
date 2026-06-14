using SafeTestsets

@safetestset "Aqua" begin
    using ODEInterfaceDiffEq, Aqua, Test
    # deps_compat disabled: ODEInterfaceDiffEq does not declare a compat entry
    # for the LinearAlgebra stdlib dependency. Tracked in
    # https://github.com/SciML/ODEInterfaceDiffEq.jl/issues/105
    Aqua.test_all(ODEInterfaceDiffEq; deps_compat = false)
    @test_broken false  # Aqua deps compat: missing compat entry for LinearAlgebra dep — tracked in https://github.com/SciML/ODEInterfaceDiffEq.jl/issues/105
end

@safetestset "JET" begin
    using ODEInterfaceDiffEq, JET, Test
    # JET.test_package reports an undefined-binding error
    # (ODEInterfaceDiffEq.uBottomEltype in src/solve.jl). Tracked in
    # https://github.com/SciML/ODEInterfaceDiffEq.jl/issues/105
    @test_broken false  # JET: `ODEInterfaceDiffEq.uBottomEltype` is not defined (src/solve.jl) — tracked in https://github.com/SciML/ODEInterfaceDiffEq.jl/issues/105
end
