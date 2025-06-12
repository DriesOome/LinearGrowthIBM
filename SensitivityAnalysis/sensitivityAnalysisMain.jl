include("./nominalParameters.jl")
include("./sensitivityAnalysis.jl")

### Analysis: local sensitivy analysis 
## step 1: create ranges for each parameter under investigation
saRanges::Dict{Symbol, Vector{Float64}} = Dict{Symbol, Vector{Float64}}()
saRanges[:startingEssentialProteinConcentration] = [1.0, 400.0]
saRanges[:muMax] = [0.5, 3.0]
saRanges[:essentialProteinProductionRate] = [0.0, 5.0]
saRanges[:essentialProteinDegradationRate] = [0.0, 1.0]
saRanges[:essentialMetaboliteProductionRate] = [0.0, 5.0]
saRanges[:essentialMetaboliteDegradationRate] = [0.0, 5.0]
saRanges[:essentialMetaboliteKm] = [5.0, 100.0]
saRanges[:essentialMetaboliteThreshold] = [5.0, 100.0]



## step 2: Create and run bioreactors 
saBioreactors::Dict{Symbol, Vector{Bioreactor}} = Dict{Symbol, Vector{Bioreactor}}()
for varName in keys(saRanges)
    bioreactors::Vector{Bioreactor} = constructBioreactors(constructParameterRange(nominalParameters, varName, saRanges[varName], 10))
    runBioreactors(bioreactors)
    saBioreactors[varName] = bioreactors
end

## step 3: Plot results
savefig(plotSensitivityAnalysis(saBioreactors), "./saAnalysis.png")
