using Measures
using Plots
include("../bioreactor.jl")

# single plots
function plotBiomass(bioreactor::Bioreactor)
    return plotBiomass!(plot(), bioreactor)
end

function plotBiomass!(parentPlot, bioreactor::Bioreactor)
    fontsize=14

    # mappings
    time = collect(0:bioreactor.parameters.agentTimeStep:bioreactor.parameters.duration)

    parentPlot = plot!(parentPlot, time, [sum(bioreactor.solution(t)[getCellIdx():3:end]) for t in time], labels="total biomass", c=:blue, xlabel="time (h)", ylabel="density/L")

    parentPlot = plot!(parentPlot, xguidefontsize=fontsize, yguidefontsize=fontsize, legendfontpointsize=fontsize, legend=:none)
    parentPlot = plot!(parentPlot, xtick=[i for i in 0:1:100])
    return parentPlot
end

function plotCellCounts(bioreactor::Bioreactor)
    return plotCellCounts!(plot(), bioreactor)
end

function plotCellCounts!(parentPlot, bioreactor::Bioreactor)
    fontsize=14

    # mappings
    time = collect(0:bioreactor.parameters.agentTimeStep:bioreactor.parameters.duration)

    parentPlot = plot!(parentPlot, time, [totalCells(bioreactor.solution(t)) for t in time], labels="total cell count", c=:blue, xlabel="time (h)", ylabel="cells/L")

    parentPlot = plot!(parentPlot, xguidefontsize=fontsize, yguidefontsize=fontsize, legendfontpointsize=fontsize, legend=:none)
    parentPlot = plot!(parentPlot, xtick=[i for i in 0:1:100])
    return parentPlot
end

function plotGrowingFraction(bioreactor::Bioreactor)
    return plotGrowingFraction!(plot(), bioreactor)
end

function plotGrowingFraction!(parentPlot, bioreactor::Bioreactor)
    fontsize=14

    # mappings
    time = collect(0:bioreactor.parameters.agentTimeStep:bioreactor.parameters.duration)

    parentPlot = plot!(parentPlot, time, [100*getFractionGrowingCells(bioreactor.solution(t), bioreactor.parameters.essentialMetaboliteThreshold) for t in time], labels="fraction growing", c=:blue, xlabel="time (h)", ylabel="% growing")

    parentPlot = plot!(parentPlot, xguidefontsize=fontsize, yguidefontsize=fontsize, legendfontpointsize=fontsize, legend=:none)
    parentPlot = plot!(parentPlot, xtick=[i for i in 0:1:100], ylims=[0,105])
    return parentPlot
end

function plotLogBiomass(bioreactor::Bioreactor)

    biomassPlot = plotBiomass(bioreactor)
    biomassPlot = plot!(biomassPlot, yaxis=(:log10, [10^0, :auto]), ytick=[10^x for x in 0:12],xtick=[i for i in 0:1:100], ylabel="log10(density/L)")
    return biomassPlot
end

function plotLogBiomass!(parentPlot, bioreactor::Bioreactor)

    parentPlot = plotBiomass!(parentPlot, bioreactor)
    parentPlot = plot!(parentPlot, yaxis=(:log10, [10^0, :auto]), ytick=[10^x for x in 0:12],xtick=[i for i in 0:1:100], ylabel="log10(density/L)")
    return parentPlot
end

function plotEssentialProteinPerCell(bioreactor::Bioreactor)
    return plotEssentialProteinPerCell!(plot(), bioreactor)
end

function plotEssentialProteinPerCell!(parentPlot, bioreactor::Bioreactor)
    fontsize=14

    # mappings
    time = collect(0:bioreactor.parameters.agentTimeStep:bioreactor.parameters.duration)

    parentPlot = plot!(parentPlot, time, [sum(bioreactor.solution(t)[getCellIdx()+1:3:end]./bioreactor.solution(t)[getCellIdx():3:end])/totalCells(bioreactor.solution(t)) for t in time], c=:blue, xlabel="time (h)", ylabel="[E]/cell")

    parentPlot = plot!(parentPlot, xguidefontsize=fontsize, yguidefontsize=fontsize, legendfontpointsize=fontsize, legend=:none)
    parentPlot = plot!(parentPlot, yaxis=(:log10, [10^0, :auto]), ytick=[10^x for x in 0:12],xtick=[i for i in 0:1:100])
    return parentPlot
end

function plotEssentialMetabolitePerCell(bioreactor::Bioreactor)
    return plotEssentialMetabolitePerCell!(plot(), bioreactor)
end

function plotEssentialMetabolitePerCell!(parentPlot, bioreactor::Bioreactor)
    fontsize=14

    # mappings
    time = collect(0:bioreactor.parameters.agentTimeStep:bioreactor.parameters.duration)

    parentPlot = plot!(parentPlot, time, [sum(bioreactor.solution(t)[getCellIdx()+2:3:end]./bioreactor.solution(t)[getCellIdx():3:end])/totalCells(bioreactor.solution(t)) for t in time], c=:blue, xlabel="time (h)", ylabel="[M]/cell")

    parentPlot = plot!(parentPlot, xguidefontsize=fontsize, yguidefontsize=fontsize, legendfontpointsize=fontsize, legend=:none)
    parentPlot = plot!(parentPlot, yaxis=(:log10, [10^0, :auto]), ytick=[10^x for x in 0:12],xtick=[i for i in 0:1:100])
    return parentPlot
end

# composite plots
function plotBioreactor(bioreactor::Bioreactor)
    l = @layout [grid(1,2){0.75h}; d e]

    biomassPlot = plotBiomass(bioreactor)
    logBiomassPlot = plotLogBiomass(bioreactor)
    essentialProteinPlot = plotEssentialProteinPerCell(bioreactor)
    essentialMetabolitePlot = plotEssentialMetabolitePerCell(bioreactor)

    return plot!(biomassPlot, logBiomassPlot, essentialProteinPlot, essentialMetabolitePlot, layout=l, size=(2560,1444).*0.7, margin=10mm)
end

function plotBioreactors(bioreactors::Vector{Bioreactor})
    fontsize=14
    biomassPlot = plot()
    essentialProteinPlot = plot()
    essentialMetabolitePlot = plot()

    for bioreactor in bioreactors
        time = collect(0:bioreactor.parameters.agentTimeStep:bioreactor.parameters.duration)
        biomassPlot = plot!(biomassPlot, time, [sum(bioreactor.solution(t)[getCellIdx():3:end]) for t in time], labels="total biomass", xlabel="time (h)", ylabel="density/L")
        essentialProteinPlot = plot!(essentialProteinPlot, time, [sum(bioreactor.solution(t)[getCellIdx()+1:3:end]./bioreactor.solution(t)[getCellIdx():3:end])/totalCells(bioreactor.solution(t)) for t in time], labels="[E]/cell", xlabel="time (h)", ylabel="E")
        essentialMetabolitePlot = plot!(essentialMetabolitePlot, time, [sum(bioreactor.solution(t)[getCellIdx()+2:3:end]./bioreactor.solution(t)[getCellIdx():3:end])/totalCells(bioreactor.solution(t)) for t in time], labels="[M]/cell", xlabel="time (h)", ylabel="M")
    end
    essentialProteinPlot = plot!(essentialProteinPlot, xguidefontsize=fontsize, yguidefontsize=fontsize, legendfontpointsize=fontsize, legend=:best)
    essentialProteinPlot = plot!(essentialProteinPlot, yaxis=(:log10, [10^0, :auto]), ytick=[10^x for x in 0:12],xtick=[i for i in 0:1:100])
    biomassPlot = plot!(biomassPlot, xguidefontsize=fontsize, yguidefontsize=fontsize, legendfontpointsize=fontsize, legend=:bottomright)
    biomassPlot = plot!(biomassPlot,yaxis=(:log10, [10^0, :auto]), ytick=[10^x for x in 0:12],xtick=[i for i in 0:1:100])
    essentialMetabolitePlot = plot!(essentialMetabolitePlot, xguidefontsize=fontsize, yguidefontsize=fontsize, legendfontpointsize=fontsize, legend=:best)
    essentialMetabolitePlot = plot!(essentialMetabolitePlot, yaxis=(:log10, [10^0, :auto]), ytick=[10^x for x in 0:12],xtick=[i for i in 0:1:100])

    l = @layout [a{0.75h}; d e]
    return plot!(biomassPlot, essentialProteinPlot, essentialMetabolitePlot, layout=l, size=(2560,1444).*0.7, margin=10mm)
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

function calculateDerivative(timepoints::Vector{Float64}, values::Vector{Float64})
    itp = interpolate((timepoints,), values, Gridded(Linear()))
    return only.(Interpolations.gradient.(Ref(itp), timepoints))
end

# histogram plots
function essentialProteinHistogram(bioreactor::Bioreactor, timepoint::Float64)
    essentialProteinConcentrations::Vector{Float64} = bioreactor.solution(timepoint)[getCellIdx()+1:3:end]
    return histogram(essentialProteinConcentrations, normalize=:probability)
end

function essentialMetaboliteHistogram(bioreactor::Bioreactor, timepoint::Float64)
    essentialMetaboliteConcentrations::Vector{Float64} = bioreactor.solution(timepoint)[getCellIdx()+2:3:end]
    return histogram(essentialMetaboliteConcentrations, normalize=:probability)
end