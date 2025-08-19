include("../bioreactorParameters.jl")
include("../bioreactor.jl")
using Interpolations
using Measures
using Plots
using DataFrames
using ProgressBars

# mutater functions
function mutateParameter(value::Float64, range::Vector{Float64})::Float64
    newValue::Float64 = rand(Uniform(range[1], range[2])) # value*exp(rand(Normal(0, 1)))
    if newValue < range[1] || newValue > range[2]
        return mutateParameter(value, range)
    else
        return newValue
    end
end

function mutateParameterSet(oldParameters::BioreactorParameters, saRanges::Dict{Symbol, Vector{Float64}})::BioreactorParameters
    newParameters::BioreactorParameters = deepcopy(oldParameters)
    for varName in keys(saRanges)
        if rand(Bernoulli(0.25))
            setproperty!(newParameters, varName, mutateParameter(getproperty(newParameters, varName), saRanges[varName]))
        end
    end
    return newParameters
end

# score functions
function calculateDerivative(timepoints::Vector{Float64}, values::Vector{Float64})
    itp = interpolate((timepoints,), values, Gridded(Linear()))
    return only.(Interpolations.gradient.(Ref(itp), timepoints))
end

function calculateSaMeasure(b::Bioreactor)::Float64
    timepoints = collect(0:b.parameters.agentTimeStep:b.parameters.duration)
    biomass = [sum(b.solution(t)[getCellIdx():3:end]) for t in timepoints]
    k = min(10, length(timepoints))
    timepoints = [mean(timepoints[(i-k):(i+k)]) for i in (1+k):(length(timepoints)-k)]
    biomass = [mean(biomass[(i-k):(i+k)]) for i in (1+k):(length(biomass)-k)]
    biomass ./= maximum(biomass)
    dBdt = calculateDerivative(timepoints, biomass)
    dbdt = calculateDerivative(timepoints, dBdt)
    return dbdt[end]/maximum(dbdt)
end

function calculateAcceptanceProbability(newMeasure::Float64, oldMeasure::Float64, sensitivityConstant::Float64)::Float64
    return min(1.0, exp(-(newMeasure^2)/(sensitivityConstant^2))/exp(-(oldMeasure^2)/(sensitivityConstant^2)))
end

# data functions 
function getSaBioreactorData(b::Bioreactor, saRanges::Dict{Symbol, Vector{Float64}})::Dict{Symbol, Float64}
    bData::Dict{Symbol, Float64} = Dict{Symbol, Float64}()
    for varName in keys(saRanges)
        bData[varName] = getproperty(b.parameters, varName)
    end
    bData[:measure] = calculateSaMeasure(b)
    return bData
end

# monte carlo 
function runMarkovChainMonteCarloSa(nominalParameters::BioreactorParameters, saRanges::Dict{Symbol, Vector{Float64}}, totalRuns::Int64, sensitivityConstant::Float64)::DataFrame 
    saData::DataFrame = DataFrame()
    # init first bioreactor 
    bioreactor::Bioreactor = Bioreactor(nominalParameters)
    simulateBioreactor(bioreactor)
    append!(saData, getSaBioreactorData(bioreactor, saRanges))

    # run iterations
    currentParameters::BioreactorParameters = deepcopy(nominalParameters)
    pb::ProgressBar = ProgressBar(total=totalRuns)
    for runId in 1:totalRuns
        proposedParameters::BioreactorParameters = mutateParameterSet(currentParameters, saRanges)
        try
            bioreactor = Bioreactor(proposedParameters)
            simulateBioreactor(bioreactor)
        catch e end

        bData::Dict{Symbol, Float64} = getSaBioreactorData(bioreactor, saRanges)
        if rand(Bernoulli(calculateAcceptanceProbability(bData[:measure], saData[end, :measure], sensitivityConstant)))
            saData = append!(saData, getSaBioreactorData(bioreactor, saRanges))
            currentParameters = proposedParameters
            # TEMPPPP
            timepoints = collect(0:5/60:bioreactor.parameters.duration)
            biomass = [sum(bioreactor.solution(t)[getCellIdx():3:end]) for t in timepoints]
            display(plot(timepoints, biomass))
        end
        update(pb)
    end
    return saData
end


# sa plots 
function plotSaParameters(saData::DataFrame)
    saParameterPlots = []
    for colName in names(saData)
        saParameterPlot = histogram(saData[:, colName], title=colName, legend=:none)
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


