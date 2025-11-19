
include("./mcmcBaySa.jl")
include("./nominalParameters.jl")
include("./mcmcBayPlotting.jl")
using CSV
### Main script for running the main mcmc baysian sensitivity analysis 
## step 1: create ranges for each parameter under investigation
saRanges::Dict{Symbol, Vector{Float64}} = Dict{Symbol, Vector{Float64}}()
saRanges[:muMax] = [0.5, 2.0]
saRanges[:essentialMetaboliteProductionRate] = [0.1, 5.0]
saRanges[:essentialMetaboliteDegradationRate] = [0.1, 5.0]
saRanges[:essentialMetaboliteKm] = [1e1, 1e3]

## step 2: run mcmc 
sigma = 0.05
saData = runMarkovChainMonteCarloSa(nominalParameters, saRanges, 1000, sigma, 4)
CSV.write("./SensitivityAnalysis/saData.csv", saData)

## step 3: plot results
saData = CSV.read("./SensitivityAnalysis/saData.csv", DataFrame)
burnin::Int64 = 100
saData = omitBurnin(saData, burnin)

combinedDiagnogstics = mcmcDiagnostics(saData, maxlag=100)

lbs::Dict{String, LaTeXString} = Dict{String, LaTeXString}()
lbs["essentialMetaboliteDegradationRate"] = L"\gamma_M"
lbs["essentialMetaboliteKm"] = L"K_M"
lbs["muMax"] = L"\mu_{max}"
lbs["essentialMetaboliteProductionRate"] = L"\alpha_M"

display(plotSaParameters(saData, lbs))
savefig(plotSaParameters(saData, lbs), "./SensitivityAnalysis/saMcmcBayAnalysis.png")

display(plotMcmcChains(saData, lbs))
savefig(plotMcmcChains(saData, lbs), "./SensitivityAnalysis/saMcmcBayTrace.png")
