mutable struct DEOptions{SType,CType}
    saveat::SType
    save_everystep::Bool
    callback::CType
end

mutable struct ODEInterfaceIntegrator{oType} <: AbstractODEIntegrator
    opts::oType
end
