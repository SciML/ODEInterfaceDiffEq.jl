using ODEInterfaceDiffEq, Aqua, JET
using Test

@testset "Aqua" begin
    Aqua.test_all(ODEInterfaceDiffEq)
end

@testset "JET" begin
    JET.test_package(ODEInterfaceDiffEq; target_defined_modules = true)
end
