using SciMLTesting, ODEInterfaceDiffEq, Test
using JET

# all_qualified_accesses_via_owners: names accessed via SciMLBase that SciMLBase does
# not own. SciMLStructures (the SciMLStructures-owned module) and recursive_bottom_eltype
# (RecursiveArrayTools-owned) are both reached through the SciMLBase reexport.
const QUALIFIED_VIA_OWNERS_IGNORE = (
    :SciMLStructures, :recursive_bottom_eltype,
)

# all_qualified_accesses_are_public: names still non-public in the registered releases.
# - SciMLBase-owned solver-extension internals (SciMLBase 3.28.1 still private).
# - recursive_bottom_eltype (RecursiveArrayTools-owned) / SciMLStructures
#   (SciMLStructures-owned), accessed via the SciMLBase reexport: non-public there.
# - DiffEqBase-owned internals (DiffEqBase 7.6.0 still private).
# - ODEInterface C-wrapper solver entry points and constants (ODEInterface 0.5.1, no
#   `public` declarations).
const QUALIFIED_ARE_PUBLIC_IGNORE = (
    :OUTPUTFCN_CALL_REASON, :OUTPUTFCN_CALL_STEP, :OUTPUTFCN_DENSE,
    :OUTPUTFCN_RET_CONTINUE, :OUTPUTFCN_RET_CONTINUE_XCHANGED, :OUTPUTFCN_WODENSE,
    :OptionsODE, :RHS_CALL_INSITU, :SciMLStructures, :Stats, :__solve,
    :_process_verbose_param, :calculate_solution_errors!, :ddeabm, :ddebdf, :dop853,
    :dopri5, :initialize_dae!, :odex, :radau, :radau5, :recursive_bottom_eltype,
    :rodas, :seulex, :solution_new_retcode,
)

run_qa(
    ODEInterfaceDiffEq;
    explicit_imports = true,
    ei_kwargs = (;
        all_qualified_accesses_via_owners = (; ignore = QUALIFIED_VIA_OWNERS_IGNORE),
        all_qualified_accesses_are_public = (; ignore = QUALIFIED_ARE_PUBLIC_IGNORE),
    ),
)
