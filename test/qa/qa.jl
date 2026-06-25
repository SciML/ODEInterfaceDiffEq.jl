using SciMLTesting, ODEInterfaceDiffEq, Test
using JET

# ExplicitImports ignore-lists. Every ignored name is a non-public name owned by a
# dependency (not by ODEInterfaceDiffEq), or a SciMLBase name reexported through
# DiffEqBase. These are the standard solver-extension surface: ODEInterfaceDiffEq
# extends/uses DiffEqBase/SciMLBase internals (the `DiffEqBase.__solve` /
# `initialize_dae!` convention shared by every SciML solver package) and the
# ODEInterface C-wrapper's non-public solver entry points. They become public as the
# base libraries declare them `public`; until then they are legitimately ignored.

# all_qualified_accesses_via_owners: names accessed through DiffEqBase that SciMLBase
# owns (DiffEqBase reexports them), plus SciMLStructures accessed through SciMLBase.
const QUALIFIED_VIA_OWNERS_IGNORE = (
    :AbstractODEAlgorithm, :AbstractODEIntegrator, :AbstractODEProblem,
    :AbstractParameterizedFunction, :SciMLStructures, :__solve, :build_solution,
    :calculate_solution_errors!, :has_analytic, :has_jac, :has_tgrad,
    :initialize_dae!, :solution_new_retcode,
)

# all_qualified_accesses_are_public: non-public names from DiffEqBase, SciMLBase,
# SciMLBase.ReturnCode, ODEInterface, and Base.Iterators (`filter`).
const QUALIFIED_ARE_PUBLIC_IGNORE = (
    :AbstractODEAlgorithm, :AbstractODEIntegrator, :AbstractODEProblem,
    :AbstractParameterizedFunction, :CallbackCache, :Default, :DtLessThanMin,
    :Failure, :InitialFailure, :MaxIters, :OUTPUTFCN_CALL_REASON,
    :OUTPUTFCN_CALL_STEP, :OUTPUTFCN_DENSE, :OUTPUTFCN_RET_CONTINUE,
    :OUTPUTFCN_RET_CONTINUE_XCHANGED, :OUTPUTFCN_WODENSE, :OptionsODE,
    :RHS_CALL_INSITU, :SciMLStructures, :Stats, :Success, :Unstable, :__solve,
    :_process_verbose_param, :alg_order, :apply_callback!, :apply_discrete_callback!,
    :build_solution, :calculate_solution_errors!, :ddeabm, :ddebdf, :dop853, :dopri5,
    :filter, :find_first_continuous_callback, :get_initial_values, :has_analytic,
    :has_jac, :has_tgrad, :initialize_dae!, :max_vector_callback_length, :odex,
    :radau, :radau5, :rodas, :seulex, :solution_new_retcode,
)

# all_explicit_imports_are_public: SciMLBase initialization-algorithm internals.
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
