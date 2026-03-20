using ODEInterfaceDiffEq, DiffEqBase
using Test

@time @testset "Explicit Imports" begin
    include("explicit_imports.jl")
end
@time @testset "Algorithms" begin
    include("algorithm_tests.jl")
end
@time @testset "Saving" begin
    include("saving_tests.jl")
end
@time @testset "Mass Matrix" begin
    include("mass_matrix_tests.jl")
end
@time @testset "Jacobian Tests" begin
    include("jac_tests.jl")
end
@time @testset "Callback Tests" begin
    include("callbacks.jl")
end
@time @testset "Initialization Tests" begin
    include("initialization_tests.jl")
end
@time @testset "MTK Initialization Tests" begin
    include("mtk_initialization_tests.jl")
end
