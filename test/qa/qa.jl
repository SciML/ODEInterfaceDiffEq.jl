using SciMLTesting, ODEInterfaceDiffEq, Test
using JET

# ExplicitImports ignore-lists. Every ignored name is a non-public name owned by a
# dependency (not by ODEInterfaceDiffEq): the standard SciML solver-extension surface
# (DiffEqBase/SciMLBase internals such as `__solve` / `initialize_dae!`) and the
# ODEInterface C-wrapper's non-public solver entry points. They become removable as the
# base libraries declare them `public`; until then they are legitimately ignored.

# all_qualified_accesses_via_owners: names accessed through DiffEqBase that SciMLBase
# owns (DiffEqBase reexports them), plus SciMLStructures accessed through SciMLBase, and
# recursive_bottom_eltype (RecursiveArrayTools-owned, accessed via SciMLBase reexport).
const QUALIFIED_VIA_OWNERS_IGNORE = (
    :AbstractODEAlgorithm, :AbstractODEIntegrator, :AbstractODEProblem,
    :AbstractParameterizedFunction, :SciMLStructures, :__solve, :build_solution,
    :calculate_solution_errors!, :has_analytic, :has_jac, :has_tgrad,
    :initialize_dae!, :recursive_bottom_eltype, :solution_new_retcode,
)

# all_qualified_accesses_are_public: names still non-public in the registered releases
# (DiffEqBase 7, SciMLBase 3.24, ODEInterface 0.5).
const QUALIFIED_ARE_PUBLIC_IGNORE = (
    :AbstractODEAlgorithm, :AbstractODEIntegrator, :AbstractODEProblem,
    :AbstractParameterizedFunction, :CallbackCache, :OUTPUTFCN_CALL_REASON,
    :OUTPUTFCN_CALL_STEP, :OUTPUTFCN_DENSE, :OUTPUTFCN_RET_CONTINUE,
    :OUTPUTFCN_RET_CONTINUE_XCHANGED, :OUTPUTFCN_WODENSE, :OptionsODE,
    :RHS_CALL_INSITU, :SciMLStructures, :Stats, :__solve, :_process_verbose_param,
    :alg_order, :apply_callback!, :apply_discrete_callback!, :build_solution,
    :calculate_solution_errors!, :ddeabm, :ddebdf, :dop853, :dopri5,
    :find_first_continuous_callback, :get_initial_values, :has_analytic, :has_jac,
    :has_tgrad, :initialize_dae!, :max_vector_callback_length, :odex, :radau,
    :radau5, :recursive_bottom_eltype, :rodas, :seulex, :solution_new_retcode,
)

# all_explicit_imports_are_public: SciMLBase initialization-algorithm internals,
# still non-public in SciMLBase 3.24.
const EXPLICIT_IMPORTS_ARE_PUBLIC_IGNORE = (
    :NoInit, :OverrideInit, :has_initialization_data,
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
