mutable struct DEOptions{SType,CType}
    saveat::SType
    save_on::Bool
    save_everystep::Bool
    callback::CType
end

mutable struct ODEInterfaceIntegrator{uType,uPrevType,oType,SType,solType,algType} <: DiffEqBase.AbstractODEIntegrator
    u::uType
    uprev::uPrevType
    t::Float64
    tprev::Float64
    opts::oType
    u_modified::Bool
    tdir::Float64
    sizeu::SType
    sol::solType
    eval_sol_fcn
    event_last_time::Int
    alg::algType
end

(integrator::ODEInterfaceIntegrator)(t) = integrator.eval_sol_fcn(t)
