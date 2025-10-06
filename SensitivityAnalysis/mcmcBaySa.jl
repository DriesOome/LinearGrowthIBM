using Distributed
addprocs(max(4-nprocs(), 0))

@everywhere include("../bioreactorParameters.jl")
@everywhere include("../bioreactor.jl")
@everywhere using Interpolations
using Measures
using Plots
@everywhere using DataFrames
using ProgressBars

# mutater functions
@everywhere function mutateParameter(value::Float64, range::Vector{Float64})::Float64
    newValue::Float64 = rand(Uniform(range[1], range[2])) # value*exp(rand(Normal(0, 1)))
    if newValue < range[1] || newValue > range[2]
        return mutateParameter(value, range)
    else
        return newValue
    end
end

@everywhere function mutateParameterSet(oldParameters::BioreactorParameters, saRanges::Dict{Symbol, Vector{Float64}})::BioreactorParameters
    newParameters::BioreactorParameters = deepcopy(oldParameters)
    for varName in keys(saRanges)
        setproperty!(newParameters, varName, mutateParameter(getproperty(newParameters, varName), saRanges[varName]))
    end
    return newParameters
end

@everywhere function enforceParameterInvariants(parameters::BioreactorParameters)
    # steady state [M] should give 
    return parameters
end

# score functions
@everywhere function calculateDerivative(timepoints::Vector{Float64}, values::Vector{Float64})
    itp = interpolate((timepoints,), values, Gridded(Linear()))
    return only.(Interpolations.gradient.(Ref(itp), timepoints))
end

@everywhere function calculateSaMeasure(b::Bioreactor)::Float64
    timepoints = collect(0:b.parameters.agentTimeStep:b.parameters.duration)
    biomass = [getTotalCellDensity(b.solution(t)) for t in timepoints]
    k = min(10, length(timepoints))
    timepoints = [mean(timepoints[(i-k):(i+k)]) for i in (1+k):(length(timepoints)-k)]
    biomass = [mean(biomass[(i-k):(i+k)]) for i in (1+k):(length(biomass)-k)]
    dBdt = calculateDerivative(timepoints, biomass)
    dbdt = calculateDerivative(timepoints, dBdt)./(b.parameters.startingCellCount*b.parameters.muMax^2)
    return dbdt[end]/maximum(dbdt)
end

@everywhere function calculateAcceptanceProbability(newMeasure::Float64, oldMeasure::Float64, sensitivityConstant::Float64)::Float64
    return min(1.0, exp((oldMeasure^2-newMeasure^2)/sensitivityConstant^2))
    #return min(1.0, exp(-(newMeasure^2)/(sensitivityConstant^2))/exp(-(oldMeasure^2)/(sensitivityConstant^2)))
end

# data functions 
@everywhere function getSaBioreactorData(b::Bioreactor, saRanges::Dict{Symbol, Vector{Float64}})::Dict{Symbol, Float64}
    bData::Dict{Symbol, Float64} = Dict{Symbol, Float64}()
    for varName in keys(saRanges)
        bData[varName] = getproperty(b.parameters, varName)
    end
    bData[:measure] = calculateSaMeasure(b)
    return bData
end


# monte carlo 
@everywhere function runMarkovChainMonteCarloSa(nominalParameters::BioreactorParameters, saRanges::Dict{Symbol, Vector{Float64}}, totalRuns::Int64, sensitivityConstant::Float64)::DataFrame 
    saData::DataFrame = DataFrame()
    # init first bioreactor 
    bioreactor::Bioreactor = Bioreactor(nominalParameters)
    simulateBioreactor(bioreactor)
    bData = getSaBioreactorData(bioreactor, saRanges)
    currentMeasure = bData[:measure]
    append!(saData, bData)

    # run iterations
    currentParameters::BioreactorParameters = deepcopy(nominalParameters)
    pb::ProgressBar = ProgressBar(total=totalRuns)
    for runId in 1:totalRuns
        proposedParameters::BioreactorParameters = mutateParameterSet(currentParameters, saRanges)
        try
            bioreactor = Bioreactor(proposedParameters)
            simulateBioreactor(bioreactor)
            bData::Dict{Symbol, Float64} = getSaBioreactorData(bioreactor, saRanges)
            if rand(Bernoulli(calculateAcceptanceProbability(bData[:measure], currentMeasure, sensitivityConstant)))
                currentParameters = proposedParameters
                currentMeasure = bData[:measure]
            end
            saData = append!(saData, getSaBioreactorData(bioreactor, saRanges))
        catch e 
            println(e)
        end
        update(pb)
    end
    return saData
end

function runMarkovChainMonteCarloSa(nominalParameters::BioreactorParameters, saRanges::Dict{Symbol, Vector{Float64}}, totalRuns::Int64, sensitivityConstant::Float64, totalChains::Int64)::DataFrame
    # run chains
    chains = pmap(1:totalChains) do chainId
        runMarkovChainMonteCarloSa(nominalParameters, saRanges, totalRuns, sensitivityConstant)
    end
    # collect chains
    saData::DataFrame = DataFrame()
    for chainId in 1:totalChains
        saSubData::DataFrame = chains[chainId]
        saSubData[!, :chainId] .= chainId
        append!(saData, saSubData) 
    end
    return saData
end

# sa plots 
function plotSaParameters(saData::DataFrame)
    saParameterPlots = []
    for colName in names(saData)
        saParameterPlot = histogram(saData[:, colName], title=colName, legend=:none, normalize=:pdf)
        push!(saParameterPlots, saParameterPlot)
    end
    return plot(saParameterPlots..., size=(2560,1444).*0.7, margin=10mm)
end

# derivative plots
function plotBiomassDerivative(bioreactor::Bioreactor)
    return plotBiomassDerivative!(plot(), bioreactor)
end

function plotBiomassDerivative!(parentPlot, bioreactor::Bioreactor)
    timepoints::Vector{Float64} = collect(0:bioreactor.parameters.agentTimeStep:bioreactor.parameters.duration)
    biomass::Vector{Float64} = [sum(bioreactor.solution(t)[getCellIdx():3:end]) for t in timepoints]
    derivative = calculateDerivative(timepoints, biomass)
    return plot!(parentPlot, timepoints, derivative)
end

function plotMcmcChains(saData::DataFrame)
    saParameterPlots = []
    for colName in names(saData)
        saParameterPlot = plot(1:size(saData, 1), saData[:, colName], title=colName, legend=:none)
        push!(saParameterPlots, saParameterPlot)
    end
    return plot(saParameterPlots..., size=(2560,1444).*0.7, margin=10mm)
end

function omitBurnin(saData::DataFrame, burnin::Int64)
    burninData = DataFrame()
    for saSubData in groupby(saData, :chainId)
        saSubData = saSubData[Not(1:burnin), :]
        saSubData[!, :iteration] = collect(1:nrow(saSubData))
        append!(burninData, saSubData)
    end
    return burninData 
end

function plotAutoCorrelation(saData::DataFrame; maxlag=100)
    saAutoCorrelationPlots = []
    for colName in setdiff(names(saData), ["measure", "iteration", "chainId"])
        saAutocorrelationPlot = plot()
        for subData in groupby(saData, :chainId)
            plot!(saAutocorrelationPlot, 0:maxlag, autocor(subData[:, colName], 0:maxlag), title=colName, legend=:none)
        end
        push!(saAutoCorrelationPlots, saAutocorrelationPlot)
    end
    return plot(saAutoCorrelationPlots..., size=(2560,1444).*0.7, margin=10mm)
end


"""
Compute autocorrelation, IAT, and ESS for each column of a DataFrame of MCMC samples.
Returns a DataFrame with one row per parameter.
"""
function chainDiagnostics(df::Union{DataFrame, SubDataFrame}; maxlag::Int=2000, window_rule=:firstneg)
    N = nrow(df)
    params = setdiff(names(saData), ["chainId", "measure", "iteration"])
    result = DataFrame(parameter=String[], iat=Float64[], ess=Float64[], ss=Int[], mu=Float64[], sigma=Float64[])
    for p in params
        x = df[:, p]
        rho = autocor(x, 0:maxlag)
        # truncate at first negative or use full sum
        idx = findfirst(<(0), rho[2:end])
        M = isnothing(idx) ? length(rho)-1 : idx
        τ = 1 + 2 * sum(rho[2:M])
        ESS = N / τ
        mu = mean(x)
        sigma = var(x)
        push!(result, (String(p), τ, ESS, length(x), mu, sigma))
    end
    return result
end

function mcmcDiagnostics(saData::DataFrame; maxlag::Int64=2000, window_rule=:firstneg)
    combinedDiagnogstics = DataFrame()
    for saSubData in groupby(saData, :chainId)
        subDiagnostics = chainDiagnostics(saSubData, maxlag=maxlag)
        subDiagnostics[!, :chainId] .= saSubData[1, :chainId]
        append!(combinedDiagnogstics, subDiagnostics)
    end
    return combine(groupby(combinedDiagnogstics, :parameter)) do sub
        IAT_combined = sum(sub.ss)/sum(sub.ss./sub.iat)
        ESS_combined = sum(sub.ess)
        n = mean(sub.ss)
        mu_g = mean(sub.mu)
        B = n*var(sub.mu)
        W = mean(sub.sigma)
        V = ((n-1)/n)*W + (1/n)*B
        (; parameter=first(sub.parameter), iat_combined = IAT_combined, combined_ess=ESS_combined, rhat=sqrt(V/W))
    end
end


