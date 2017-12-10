mutable struct DEOptions{SType,CType}
    saveat::SType
    save_everystep::Bool
    callback::CType
end

mutable struct ODEInterfaceIntegrator{uType,oType,SType,solType} <: AbstractODEIntegrator
    u::uType
    t::Float64
    tprev::Float64
    opts::oType
    u_modified::Bool
    tdir::Float64
    sizeu::SType
    sol::solType
end
