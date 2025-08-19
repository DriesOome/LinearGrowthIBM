include("./mcmcBaySa.jl")
include("./nominalParameters.jl")

### Main script for running the main mcmc baysian sensitivity analysis 
#

### Analysis: local sensitivy analysis 
## step 1: create ranges for each parameter under investigation
saRanges::Dict{Symbol, Vector{Float64}} = Dict{Symbol, Vector{Float64}}()
saRanges[:startingEssentialProteinConcentration] = [50.0, 400.0]
saRanges[:muMax] = [0.5, 2.0]
saRanges[:essentialMetaboliteProductionRate] = [0.0, 5.0]
saRanges[:essentialMetaboliteDegradationRate] = [0.0, 5.0]
saRanges[:essentialMetaboliteKm] = [5.0, 200.0]
# saRanges[:essentialMetaboliteThreshold] = [5.0, 100.0]

## step 2: run mcmc 
saData::DataFrame = runMarkovChainMonteCarloSa(nominalParameters, saRanges, 1000, 0.2)

## step 3: plot results 
plotSaParameters(saData)

