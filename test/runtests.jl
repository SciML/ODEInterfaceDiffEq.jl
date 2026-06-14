using ODEInterfaceDiffEq, DiffEqBase
using SafeTestsets, Test

const GROUP = get(ENV, "GROUP", "All")

if GROUP == "All" || GROUP == "Core"
    @time @safetestset "Explicit Imports" begin
        include("explicit_imports.jl")
    end
    @time @safetestset "Algorithms" begin
        include("algorithm_tests.jl")
    end
    @time @safetestset "Saving" begin
        include("saving_tests.jl")
    end
    @time @safetestset "Mass Matrix" begin
        include("mass_matrix_tests.jl")
    end
    @time @safetestset "Jacobian Tests" begin
        include("jac_tests.jl")
    end
    @time @safetestset "Callback Tests" begin
        include("callbacks.jl")
    end
    @time @safetestset "Initialization Tests" begin
        include("initialization_tests.jl")
    end
    @time @safetestset "MTK Initialization Tests" begin
        include("mtk_initialization_tests.jl")
    end
end

if GROUP == "QA"
    using Pkg
    Pkg.activate(joinpath(@__DIR__, "qa"))
    Pkg.instantiate()
    @time @testset "Quality Assurance" begin
        include(joinpath("qa", "qa.jl"))
    end
end
