include("../bioreactor.jl")
using Plots
using Measures
using Interpolations
using LaTeXStrings
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
        simulateBioreactor(bioreactor)
    end
end

function plotSensitivityAnalysis(saBioreactors::Dict{Symbol, Vector{Bioreactor}}, labels::Dict{Symbol, LaTeXString})
    l = @layout [a b c d; e f g h]
    bPlots = []
    for varName in keys(saBioreactors)
        bPlot = plotSaBioreactors(saBioreactors[varName], varName, labels[varName])
        push!(bPlots, bPlot)
    end
    return plot(bPlots..., layout=l, size=(1600,1000), margin=6mm)
end

function plotSaBioreactors(bioreactors::Vector{Bioreactor}, varName::Symbol, label::LaTeXString)
    biomassPlot = plot()
    for bioreactor in bioreactors 
        biomassPlot = plotSaBiomass!(biomassPlot, bioreactor, varName)
        biomassPlot = plot!(biomassPlot, title=label)
    end
    return biomassPlot
end

function plotSaBiomass!(parentPlot, bioreactor::Bioreactor, varName::Symbol)
    fontsize=14
    time = collect(0:bioreactor.parameters.agentTimeStep:bioreactor.parameters.duration)

    parentPlot = plot!(parentPlot, time, [getTotalCellDensity(bioreactor.solution(t)) for t in time], labels=string(round(getproperty(bioreactor.parameters, varName), digits=1)))
    parentPlot = plot!(parentPlot, xlabel="time (h)", ylabel="density")
    parentPlot = plot!(parentPlot, xguidefontsize=fontsize, yguidefontsize=fontsize, legendfontpointsize=7, legend=:best)
    parentPlot = plot!(parentPlot, xtick=[i for i in 0:1:100])
    return parentPlot
end