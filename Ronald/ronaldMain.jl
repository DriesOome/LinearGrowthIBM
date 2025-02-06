include("../bioreactor.jl")
include("../Plotting/plotting.jl")
include("./ronaldParameters.jl")
include("../dataExtractor.jl")
include("./ronaldPlots.jl")

dataTimePoints = collect(0:0.5:ronaldParameters.duration)

#= symmetric growth with dilution
b = Bioreactor(ronaldParameters)
simulateBioreactor(b, duration)
scenarioName::String = "exp_constraint"
saveBioreactorRun(b, dataTimePoints, "./Ronald/Data/"*scenarioName)
biomassPlotExpConstr = plotBiomass(b)
cellCountPlotExpConstr = plotCellCounts(b)
growingFractionPlotExpConstr = plotGrowingFraction(b)
savefig(plotBiomass(b), "./Ronald/Figures/"*scenarioName*"_biomassPlot.png")
savefig(plotCellCounts(b), "./Ronald/Figures/"*scenarioName*"_cellCountPlot.png")
savefig(plotGrowingFraction(b), "./Ronald/Figures/"*scenarioName*"_growingFractionPlot.png")
=#
# assymetric growth
parametersLinear = deepcopy(ronaldParameters)
parametersLinear.divisionSymmetry = 0.0
bLinear = Bioreactor(parametersLinear)
simulateBioreactor(bLinear, parametersLinear.duration)
scenarioName = "linear"
saveBioreactorRun(bLinear, dataTimePoints, "./Ronald/Data/"*scenarioName)
biomassPlotLinear = plotBiomass(bLinear)
cellCountPlotLinear = plotCellCounts(bLinear)
growingFractionPlotLinear = plotGrowingFraction(bLinear)
ronaldPlotLinear = plotRonald(bLinear)

savefig(plotBiomass(bLinear), "./Ronald/Figures/"*scenarioName*"_biomassPlot.png")
savefig(plotCellCounts(bLinear), "./Ronald/Figures/"*scenarioName*"_cellCountPlot.png")
savefig(plotGrowingFraction(bLinear), "./Ronald/Figures/"*scenarioName*"_growingFractionPlot.png")

# symmetric unlimited growth
parametersExpUnConstr = deepcopy(ronaldParameters)
parametersExpUnConstr.divisionSymmetry = 0.5
parametersExpUnConstr.startingEssentialProteinConcentration = 10e6
parametersExpUnConstr.duration = 3.0
bExpUnConstr = Bioreactor(parametersExpUnConstr)
simulateBioreactor(bExpUnConstr, parametersExpUnConstr.duration)
scenarioName = "exp_unconstraint"
saveBioreactorRun(bExpUnConstr, dataTimePoints, "./Ronald/Data/"*scenarioName)
biomassPlotExpUnConstr = plotBiomass(bExpUnConstr)
cellCountPlotExpUnConstr = plotCellCounts(bExpUnConstr)
growingFractionPlotExpUnConstr = plotGrowingFraction(bExpUnConstr)
ronaldPlotExpUnConstr = plotRonald(bExpUnConstr)
savefig(plotBiomass(bExpUnConstr), "./Ronald/Figures/"*scenarioName*"_biomassPlot.png")
savefig(plotCellCounts(bExpUnConstr), "./Ronald/Figures/"*scenarioName*"_cellCountPlot.png")
savefig(plotGrowingFraction(bExpUnConstr), "./Ronald/Figures/"*scenarioName*"_growingFractionPlot.png")

savefig(plotRonald(bExpUnConstr, bLinear), "./Ronald/Figures/globalPlot.png")



