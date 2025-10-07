using Plots
using Measures
using LaTeXStrings
# layout 
figSize = (1600, 1000)
xlabelsize = 20

# sa plots 
function plotSaParameters(saData::DataFrame, labels::Dict{String, LaTeXString})
    saParameterPlots = []
    for colName in setdiff(names(saData), ["measure", "chainId", "iteration"])
        saParameterPlot = histogram(saData[:, colName], xlabel=labels[colName], legend=:none, normalize=:pdf, labelfontsize=xlabelsize)
        push!(saParameterPlots, saParameterPlot)
    end
    pLeft = plot(saParameterPlots..., size=figSize, margin=6mm, layout=(2,2))
    pRight = histogram(saData[:, :measure], xlabel=L"s_\theta", legend=:none, normalize=:pdf, labelfontsize=xlabelsize)
    return plot(pLeft, pRight, layout=(1,2))
end

function plotSaScatter(saData::DataFrame)
    saParameterPlots = []
    parameterNames = setdiff(names(saData), ["measure", "iteration", "chainId"])
    for colName1 in parameterNames
        for colName2 in parameterNames
            zcolor = log.(saData[:, :measure].^2)
            saParameterPlot = scatter(saData[:, colName1], saData[:, colName2], cmap=:viridis, zcolor=zcolor, markersize=2, markerstrokewidth=0.0)
            push!(saParameterPlots, saParameterPlot)
        end
    end
    return plot(saParameterPlots..., size=figSize)
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

function plotMcmcChains(saData::DataFrame, labels::Dict{String, LaTeXString})
    saParameterPlots = []
    for colName in setdiff(names(saData), ["chainId", "iteration", "measure"])
        saParameterPlot = plot()
        for saSubData in groupby(saData, :chainId)
            plot!(saParameterPlot, 1:size(saSubData, 1), saSubData[:, colName],  ylabel=labels[colName],
                labelfontsize=xlabelsize, label="chain "*string(saSubData[1, :chainId]), linealpha=0.75, legend=:topright)
        end
        push!(saParameterPlots, saParameterPlot)
    end
    return plot(saParameterPlots..., size=figSize, margin=10mm, layout=(2,2))
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
