using ExplicitImports
using ODEInterfaceDiffEq
using Test

@testset "ExplicitImports" begin
    @test check_no_implicit_imports(ODEInterfaceDiffEq) === nothing
    @test check_no_stale_explicit_imports(ODEInterfaceDiffEq) === nothing
end
