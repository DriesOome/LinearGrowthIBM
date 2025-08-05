include("./mcmcBaySa.jl")
include("./nominalParameters.jl")

### Main script for running the main mcmc baysian sensitivity analysis 
#

### Analysis: local sensitivy analysis 
## step 1: create ranges for each parameter under investigation
saRanges::Dict{Symbol, Vector{Float64}} = Dict{Symbol, Vector{Float64}}()
saRanges[:startingEssentialProteinConcentration] = [1.0, 400.0]
saRanges[:muMax] = [0.5, 3.0]
saRanges[:essentialMetaboliteProductionRate] = [0.0, 5.0]
saRanges[:essentialMetaboliteDegradationRate] = [0.0, 5.0]
saRanges[:essentialMetaboliteKm] = [5.0, 100.0]
# saRanges[:essentialMetaboliteThreshold] = [5.0, 100.0]

## step 2: run mcmc 
saData::DataFrame = runMarkovChainMonteCarloSa(nominalParameters, saRanges, 1000, 0.2)

## step 3: plot results 
plotSaParameters(saData)

#= temp 
timepoints = collect(0:1:b.parameters.duration)
biomass = [sum(b.solution(t)[getCellIdx():3:end]) for t in timepoints]
dBdt = calculateDerivative(timepoints, biomass)
dbdt = calculateDerivative(timepoints, dBdt)
plot(timepoints, dBdt)
plot(timepoints, dbdt)
=#