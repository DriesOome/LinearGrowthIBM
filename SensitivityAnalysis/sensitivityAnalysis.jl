include("../bioreactor.jl")
using Plots
using Measures
using Interpolations
function constructParameterRange(nominalParameters::BioreactorParameters, varName::Symbol, range::Vector{Float64}, nSets::Int64)::Vector{BioreactorParameters}
    parameterVector::Vector{BioreactorParameters} = []
    varValues::Vector{Float64} = LinRange(range[1], range[2], nSets)
    for setId in 1:nSets
        newSet::BioreactorParameters = deepcopy(nominalParameters)
        setproperty!(newSet, varName, varValues[setId])
        push!(parameterVector, newSet)
    end
    return parameterVector
end

function constructBioreactors(parameterVector::Vector{BioreactorParameters})::Vector{Bioreactor}
    bioreactors::Vector{Bioreactor} = []
    for pSet in parameterVector
        push!(bioreactors, Bioreactor(pSet))
    end
    return bioreactors
end

function runBioreactors(bioreactors::Vector{Bioreactor})
    for bioreactor in bioreactors
        simulateBioreactor(bioreactor, bioreactor.parameters.duration)
    end
end

function plotSensitivityAnalysis(saBioreactors::Dict{Symbol, Vector{Bioreactor}})
    l = (Int64(ceil(length(saBioreactors)/4)), min(4, length(saBioreactors)))
    bPlots = []
    for varName in keys(saBioreactors)
        bPlot = plotSaBioreactors(saBioreactors[varName], varName)
        push!(bPlots, bPlot)
    end
    return plot(bPlots..., layout=l, size=(2560,1444).*0.7, margin=10mm)
end

function plotSaBioreactors(bioreactors::Vector{Bioreactor}, varName::Symbol)
    biomassPlot = plot()
    for bioreactor in bioreactors 
        biomassPlot = plotSaBiomass!(biomassPlot, bioreactor, varName)
        biomassPlot = plot!(biomassPlot, title=string(varName))
    end
    return biomassPlot
end

function plotSaBiomass!(parentPlot, bioreactor::Bioreactor, varName::Symbol)
    fontsize=14
    time = collect(0:bioreactor.parameters.agentTimeStep:bioreactor.parameters.duration)

    parentPlot = plot!(parentPlot, time, [sum(bioreactor.solution(t)[getCellIdx():3:end]) for t in time], labels=string(round(getproperty(bioreactor.parameters, varName), digits=1)))
    parentPlot = plot!(parentPlot, xlabel="time (h)", ylabel="density/L")
    parentPlot = plot!(parentPlot, xguidefontsize=fontsize, yguidefontsize=fontsize, legendfontpointsize=7, legend=:best)
    parentPlot = plot!(parentPlot, xtick=[i for i in 0:1:100])
    return parentPlot
end