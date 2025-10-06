include("./mcmcBaySa.jl")
include("./nominalParameters.jl")
using CSV
### Main script for running the main mcmc baysian sensitivity analysis 
#

### Analysis: local sensitivy analysis 
## step 1: create ranges for each parameter under investigation
saRanges::Dict{Symbol, Vector{Float64}} = Dict{Symbol, Vector{Float64}}()
#saRanges[:startingEssentialProteinConcentration] = [50.0, 400.0]
saRanges[:muMax] = [0.5, 2.0]
saRanges[:essentialMetaboliteProductionRate] = [0.1, 5.0]
saRanges[:essentialMetaboliteDegradationRate] = [0.1, 5.0]
saRanges[:essentialMetaboliteKm] = [1e1, 1e3]
# saRanges[:essentialMetaboliteThreshold] = [5.0, 100.0]

## step 2: run mcmc 
sigma = 0.01
saData = runMarkovChainMonteCarloSa(nominalParameters, saRanges, 1000, sigma, 4)
CSV.write("./SensitivityAnalysis/saData.csv", saData)

## step 3: plot results
saData = CSV.read("./SensitivityAnalysis/saData.csv", DataFrame)
burnin::Int64 = 100 
saData = omitBurnin(saData, burnin)

combinedDiagnogstics = mcmcDiagnostics(saData, maxlag=100)


display(plotSaParameters(saData))
savefig(plotSaParameters(saData), "./SensitivityAnalysis/saMcmcBayAnalysis.png")

display(plotMcmcChains(saData))
savefig(plotMcmcChains(saData), "./SensitivityAnalysis/saMcmcBayTrace.png")
