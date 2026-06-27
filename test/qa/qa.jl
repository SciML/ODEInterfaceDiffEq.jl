using SciMLTesting, ODEInterfaceDiffEq, Test
using JET

# all_qualified_accesses_via_owners: names accessed through DiffEqBase that SciMLBase
# owns (DiffEqBase reexports them), plus SciMLStructures accessed through SciMLBase, and
# recursive_bottom_eltype (RecursiveArrayTools-owned, accessed via SciMLBase reexport).
const QUALIFIED_VIA_OWNERS_IGNORE = (
    :AbstractODEAlgorithm, :AbstractODEIntegrator, :AbstractODEProblem,
    :AbstractParameterizedFunction, :SciMLStructures, :__solve, :build_solution,
    :calculate_solution_errors!, :has_analytic, :has_jac, :has_tgrad,
    :initialize_dae!, :recursive_bottom_eltype, :solution_new_retcode,
)

# all_qualified_accesses_are_public: names still non-public in the registered releases.
# - SciMLBase-owned solver-extension internals (SciMLBase 3.27.0 still private).
# - recursive_bottom_eltype (RecursiveArrayTools-owned) / SciMLStructures
#   (SciMLStructures-owned), accessed via the SciMLBase reexport: non-public there.
# - DiffEqBase-owned internals (DiffEqBase 7.6.0 still private).
# - ODEInterface C-wrapper solver entry points and constants (ODEInterface 0.5.1, no
#   `public` declarations).
const QUALIFIED_ARE_PUBLIC_IGNORE = (
    :AbstractODEIntegrator, :AbstractParameterizedFunction, :OUTPUTFCN_CALL_REASON,
    :OUTPUTFCN_CALL_STEP, :OUTPUTFCN_DENSE, :OUTPUTFCN_RET_CONTINUE,
    :OUTPUTFCN_RET_CONTINUE_XCHANGED, :OUTPUTFCN_WODENSE, :OptionsODE, :RHS_CALL_INSITU,
    :SciMLStructures, :Stats, :__solve, :_process_verbose_param,
    :calculate_solution_errors!, :ddeabm, :ddebdf, :dop853, :dopri5, :has_analytic,
    :has_tgrad, :initialize_dae!, :odex, :radau, :radau5, :recursive_bottom_eltype,
    :rodas, :seulex, :solution_new_retcode,
)

# all_explicit_imports_are_public: SciMLBase initialization-system internal, still
# non-public in SciMLBase 3.27.0.
const EXPLICIT_IMPORTS_ARE_PUBLIC_IGNORE = (
    :has_initialization_data,
)

run_qa(
    ODEInterfaceDiffEq;
    explicit_imports = true,
    ei_kwargs = (;
        all_qualified_accesses_via_owners = (; ignore = QUALIFIED_VIA_OWNERS_IGNORE),
        all_qualified_accesses_are_public = (; ignore = QUALIFIED_ARE_PUBLIC_IGNORE),
        all_explicit_imports_are_public = (; ignore = EXPLICIT_IMPORTS_ARE_PUBLIC_IGNORE),
    ),
)
