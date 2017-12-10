function solve{uType,tType,isinplace,AlgType<:ODEInterfaceAlgorithm}(
    prob::AbstractODEProblem{uType,tType,isinplace},
    alg::AlgType,
    timeseries=[],ts=[],ks=[];
    verbose=true,
    save_start=true,
    timeseries_errors=true,
    callback=nothing,kwargs...)

    isstiff = !(typeof(alg) <: Union{dopri5,dop853,odex,ddeabm})
    if verbose
        warned = !isempty(kwargs) && check_keywords(alg, kwargs, warnlist)
        if !(typeof(prob.f) <: AbstractParameterizedFunction) && isstiff
            if has_tgrad(prob.f)
                warn("Explicit t-gradient given to this stiff solver is ignored.")
                warned = true
            end
        end
        warned && warn_compat()
    end

    callbacks_internal = CallbackSet(callback,prob.callback)

    tspan = prob.tspan

    o = KW(kwargs)

    u0 = prob.u0

    if typeof(u0) <: Number
        u = [u0]
    else
        u = deepcopy(u0)
    end

    sizeu = size(u)

    if !isinplace && typeof(u)<:AbstractArray
        f! = (t,u,du) -> (du[:] = vec(prob.f(t,reshape(u,sizeu))); nothing)
    elseif !(typeof(u)<:Vector{Float64})
        f! = (t,u,du) -> (prob.f(t,reshape(u,sizeu),reshape(du,sizeu));
                          u = vec(u); du=vec(du); nothing)
    else
        f! = prob.f
    end

    o[:RHS_CALLMODE] = ODEInterface.RHS_CALL_INSITU

    _timeseries = Vector{uType}(0)
    ts = Vector{tType}(0)

    function outputfcn(reason::ODEInterface.OUTPUTFCN_CALL_REASON,
          tprev::Float64, t::Float64, u::Vector{Float64},
          eval_sol_fcn, extra_data::Dict)

      if reason == ODEInterface.OUTPUTFCN_CALL_STEP
        push!(ts,t); save_value!(_timeseries,u,uType,sizeu)
      end
      #=
      if orig_outputmode == OUTPUTFCN_NEVER || orig_outputfcn == nothing
        return OUTPUTFCN_RET_CONTINUE
      else
        return orig_outputfcn(reason,told,tnew,x,eval_sol_fcn,extra_data)
      end
      =#
      return ODEInterface.OUTPUTFCN_RET_CONTINUE
    end

    o[:OUTPUTFCN] = outputfcn
    if !(typeof(callbacks_internal.continuous_callbacks)<:Tuple{})
        o[:OUTPUTMODE] = ODEInterface.OUTPUTFCN_DENSE
    else
        o[:OUTPUTMODE] = ODEInterface.OUTPUTFCN_WODENSE
    end

    dict = buildOptions(o,
                        ODEINTERFACE_OPTION_LIST,
                        ODEINTERFACE_ALIASES,
                        ODEINTERFACE_ALIASES_REVERSED)
    if prob.mass_matrix != I
        if typeof(prob.mass_matrix) <: Matrix && isstiff
            dict[:MASSMATRIX] = prob.mass_matrix
        elseif !isstiff
            error("This solver does not support mass matrices")
        else
            error("This solver must use full or banded mass matrices.")
        end
    end
    if has_jac(prob.f)
        dict[:JACOBIMATRIX] = (t,u,J) -> prob.f(Val{:jac},t,u,J)
    end

    # Convert to the strings
    opts = ODEInterface.OptionsODE([Pair(ODEINTERFACE_STRINGS[k],v) for (k,v) in dict]...)


    if typeof(alg) <: dopri5
        _, _, retcode, stats =
            ODEInterface.dopri5(f!, tspan[1], tspan[2], vec(u), opts)
    elseif typeof(alg) <: dop853
        _, _, retcode, stats =
            ODEInterface.dop853(f!, tspan[1], tspan[2], vec(u), opts)
    elseif typeof(alg) <: odex
        _, _, retcode, stats =
            ODEInterface.odex(f!, tspan[1], tspan[2], vec(u), opts)
    elseif typeof(alg) <: seulex
        _, _, retcode, stats =
            ODEInterface.seulex(f!, tspan[1], tspan[2], vec(u), opts)
    elseif typeof(alg) <: radau
        _, _, retcode, stats =
            ODEInterface.radau(f!, tspan[1], tspan[2], vec(u), opts)
    elseif typeof(alg) <: radau5
        _, _, retcode, stats =
            ODEInterface.radau5(f!, tspan[1], tspan[2], vec(u), opts)
    elseif typeof(alg) <: rodas
        _, _, retcode, stats =
            ODEInterface.rodas(f!, tspan[1], tspan[2], vec(u), opts)
    elseif typeof(alg) <: ddeabm
        _, _, retcode, stats =
            ODEInterface.ddeabm(f!, tspan[1], tspan[2], vec(u), opts)
    elseif typeof(alg) <: ddebdf
        _, _, retcode, stats =
            ODEInterface.ddebdf(f!, tspan[1], tspan[2], vec(u), opts)
    end

    if retcode < 0
        if retcode == -1
            verbose && warn("Input is not consistent.")
            return_retcode = :Failure
        elseif retcode == -2
            verbose && warn("Interrupted. Larger maxiters is needed.")
            return_retcode = :MaxIters
        elseif retcode == -3
            verbose && warn("Step size went too small.")
            return_retcode = :DtLessThanMin
        elseif retcode == -4
            verbose && warn("Interrupted. Problem is probably stiff.")
            return_retcode = :Unstable
        end
    else
        return_retcode = :Success
    end

    build_solution(prob,  alg, ts, _timeseries,
                   timeseries_errors = timeseries_errors,
                   retcode = return_retcode)
end

function save_value!(_timeseries,u,::Type{T},sizeu) where T<:Number
    push!(_timeseries,first(u))
end

function save_value!(_timeseries,u,::Type{T},sizeu) where T<:Vector
    push!(_timeseries,u)
end

function save_value!(_timeseries,u,::Type{T},sizeu) where T<:Array
    push!(_timeseries,reshape(u,sizeu))
end

function buildOptions(o, optionlist, aliases, aliases_reversed)
    dict1 = Dict{Symbol,Any}([Pair(k,o[k]) for k in (keys(o) ∩ optionlist)])
    dict2 = Dict([Pair(aliases_reversed[k],o[k]) for k in (keys(o) ∩ values(aliases))])
    merge(dict1,dict2)
end

const ODEINTERFACE_OPTION_LIST =
    Set([:RTOL, :ATOL, :OUTPUTFCN, :OUTPUTMODE, :MAXSTEPS, :STEST, :EPS, :RHO, :SSMINSEL,
         :SSMAXSEL, :SSBETA, :MAXSS, :INITIALSS, :MAXEXCOLUMN, :STEPSIZESEQUENCE,
         :MAXSTABCHECKS, :MAXSTABCHECKLINE, :DENSEOUTPUTWOEE, :INTERPOLDEGRE,
         :SSREDUCTION, :SSSELECTPAR1, :SSSELECTPAR2, :ORDERDECFRAC, :ORDERINCFRAC,
         :OPT_RHO, :OPT_RHO2, :RHSAUTONOMOUS, :M1, :M2, :LAMBDADENSE, :TRANSJTOH,
         :STEPSIZESEQUENCE, :JACRECOMPFACTOR, :MASSMATRIX, :JACOBIMATRIX, :JACOBIBANDSSTRUCT,
         :WORKFORRHS, :WORKFORJAC, :WORKFORDEC, :WORKFORSOL,
         :MAXNEWTONITER, :NEWTONSTARTZERO, :NEWTONSTOPCRIT, :DIMFIND1VAR,
         :MAXSTAGES, :MINSTAGES, :INITSTAGES, :STEPSIZESTRATEGY,
         :FREEZESSLEFT, :FREEZESSRIGHT, :ORDERDECFACTOR,
         :ORDERINCFACTOR, :ORDERDECCSTEPFAC1, :ORDERDECSTEPFAC2, :RHS_CALLMODE
         ])

const ODEINTERFACE_ALIASES =
    Dict{Symbol,Symbol}(:RTOL=>:reltol,
                        :ATOL=>:abstol,
                        :MAXSTEPS=> :maxiters,
                        :MAXSS=>:dtmax,
                        :INITIALSS=>:dt,
                        #:SSMINSEL=>:qmin,
                        :SSBETA=>:beta2,
                        :SSMAXSEL=>:qmax)

const ODEINTERFACE_ALIASES_REVERSED =
    Dict{Symbol,Symbol}([(v,k) for (k,v) in ODEINTERFACE_ALIASES])

const ODEINTERFACE_STRINGS = Dict{Symbol,String}(
  :LOGIO            => "logio",
  :LOGLEVEL         => "loglevel",
  :RHS_CALLMODE     => "RightHandSideCallMode",

  :RTOL             => "RelTol",
  :ATOL             => "AbsTol",
  :MAXSTEPS         => "MaxNumberOfSteps",
  :EPS              => "eps",

  :OUTPUTFCN        => "OutputFcn",
  :OUTPUTMODE       => "OutputFcnMode",

  :STEST            => "StiffTestAfterStep",
  :RHO              => "rho",
  :SSMINSEL         => "StepSizeMinSelection",
  :SSMAXSEL         => "StepSizeMaxSelection",
  :SSBETA           => "StepSizeBeta",
  :MAXSS            => "MaxStep",
  :INITIALSS        => "InitialStep",


  :MAXEXCOLUMN      => "MaxExtrapolationColumn",
  :MAXSTABCHECKS    => "MaxNumberOfStabilityChecks",
  :MAXSTABCHECKLINE => "MaxLineForStabilityCheck",
  :INTERPOLDEGREE   => "DegreeOfInterpolation",
  :ORDERDECFRAC     => "OrderDecreaseFraction",
  :ORDERINCFRAC     => "OrderIncreaseFraction",
  :STEPSIZESEQUENCE => "StepSizeSequence",
  :SSREDUCTION      => "StepSizeReduction",
  :SSSELECTPAR1     => "StepSizeSelectionParam1",
  :SSSELECTPAR2     => "StepSizeSelectionParam2",
  :RHO2             => "rho2",
  :DENSEOUTPUTWOEE  => "DeactivateErrorEstInDenseOutput",

  :TRANSJTOH        => "TransfromJACtoHess",
  :MAXNEWTONITER    => "MaxNewtonIterations",
  :NEWTONSTARTZERO  => "StartNewtonWithZeros",
  :DIMOFIND1VAR     => "DimensionOfIndex1Vars",
  :DIMOFIND2VAR     => "DimensionOfIndex2Vars",
  :DIMOFIND3VAR     => "DimensionOfIndex3Vars",
  :STEPSIZESTRATEGY => "StepSizeStrategy",
  :M1               => "M1",
  :M2               => "M2",
  :JACRECOMPFACTOR  => "RecomputeJACFactor",
  :NEWTONSTOPCRIT   => "NewtonStopCriterion",
  :FREEZESSLEFT     => "FreezeStepSizeLeftBound",
  :FREEZESSRIGHT    => "FreezeStepSizeRightBound",
  :MASSMATRIX       => "MassMatrix",
  :JACOBIMATRIX     => "JacobiMatrix",
  :JACOBIBANDSTRUCT => "JacobiBandStructure",

  :MAXSTAGES        => "MaximalNumberOfStages",
  :MINSTAGES        => "MinimalNumberOfStages",
  :INITSTAGES       => "InitialNumberOfStages",
  :ORDERINCFACTOR   => "OrderIncreaseFactor",
  :ORDERDECFACTOR   => "OrderDecreaseFactor",
  :ORDERDECSTEPFAC1 => "OrderDecreaseStepFactor1",
  :ORDERDECSTEPFAC2 => "OrderDecreaseStepFactor2",

  :RHSAUTONOMOUS    => "AutonomousRHS",
  :LAMBDADENSE      => "LambdaForDenseOutput",
  :WORKFORRHS       => "WorkForRightHandSide",
  :WORKFORJAC       => "WorkForJacobimatrix",
  :WORKFORDEC       => "WorkForLuDecomposition",
  :WORKFORSOL       => "WorkForSubstitution",

  :BVPCLASS         => "BoundaryValueProblemClass",
  :SOLMETHOD        => "SolutionMethod",
  :IVPOPT           => "OptionsForIVPsolver")
