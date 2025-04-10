include("../bioreactor.jl")
include("../Plotting/plotting.jl")
include("./manuscriptParameters.jl")
include("../dataExtractor.jl")
include("./manuscriptPlots.jl")
using LaTeXStrings
dataTimePoints = collect(0:0.5:manuscriptParameters.duration)

#= symmetric growth with dilution
b = Bioreactor(manuscriptParameters)
simulateBioreactor(b, duration)
scenarioName::String = "exp_constraint"
saveBioreactorRun(b, dataTimePoints, "./ManuscriptExperiments/Data/"*scenarioName)
biomassPlotExpConstr = plotBiomass(b)
cellCountPlotExpConstr = plotCellCounts(b)
growingFractionPlotExpConstr = plotGrowingFraction(b)
savefig(plotBiomass(b), "./ManuscriptExperiments/Figures/"*scenarioName*"_biomassPlot.png")
savefig(plotCellCounts(b), "./ManuscriptExperiments/Figures/"*scenarioName*"_cellCountPlot.png")
savefig(plotGrowingFraction(b), "./ManuscriptExperiments/Figures/"*scenarioName*"_growingFractionPlot.png")
=#
# assymetric growth
parametersLinear = deepcopy(manuscriptParameters)
parametersLinear.divisionSymmetry = 0.0
bLinear = Bioreactor(parametersLinear)
simulateBioreactor(bLinear, parametersLinear.duration)
scenarioName = "linear"
saveBioreactorRun(bLinear, dataTimePoints, "./ManuscriptExperiments/Data/"*scenarioName)
biomassPlotLinear = plotBiomass(bLinear)
cellCountPlotLinear = plotCellCounts(bLinear)
growingFractionPlotLinear = plotGrowingFraction(bLinear)

savefig(plotBiomass(bLinear), "./ManuscriptExperiments/Figures/"*scenarioName*"_biomassPlot.png")
savefig(plotCellCounts(bLinear), "./ManuscriptExperiments/Figures/"*scenarioName*"_cellCountPlot.png")
savefig(plotGrowingFraction(bLinear), "./ManuscriptExperiments/Figures/"*scenarioName*"_growingFractionPlot.png")

# symmetric unlimited growth
parametersExpUnConstr = deepcopy(manuscriptParameters)
parametersExpUnConstr.divisionSymmetry = 0.5
parametersExpUnConstr.essentialProteinProductionRate = 20000.0
parametersExpUnConstr.duration = 3.0
bExpUnConstr = Bioreactor(parametersExpUnConstr)
simulateBioreactor(bExpUnConstr, parametersExpUnConstr.duration)
scenarioName = "exp_unconstraint"
saveBioreactorRun(bExpUnConstr, dataTimePoints, "./ManuscriptExperiments/Data/"*scenarioName)
biomassPlotExpUnConstr = plotBiomass(bExpUnConstr)
cellCountPlotExpUnConstr = plotCellCounts(bExpUnConstr)
growingFractionPlotExpUnConstr = plotGrowingFraction(bExpUnConstr)
savefig(plotBiomass(bExpUnConstr), "./ManuscriptExperiments/Figures/"*scenarioName*"_biomassPlot.png")
savefig(plotCellCounts(bExpUnConstr), "./ManuscriptExperiments/Figures/"*scenarioName*"_cellCountPlot.png")
savefig(plotGrowingFraction(bExpUnConstr), "./ManuscriptExperiments/Figures/"*scenarioName*"_growingFractionPlot.png")

display(plotManuscriptFigure(bExpUnConstr, bLinear))
savefig(plotManuscriptFigure(bExpUnConstr, bLinear), "./ManuscriptExperiments/Figures/globalPlot.png")



