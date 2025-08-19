using Plots
include("./bioreactor.jl")

# initial conditions
startingCellCount::Int64 = 500 # cells
startingVolume::Float64 = 400*10^-11 # L
startingEssentialProteinConcentration::Float64 = 200 # molecules

# growth kinetics
muMax::Float64 = log(2) # divisions/h
carryingCapacity::Float64 = 1e8

# essentialProtein kinetics
divisionSymmetry::Float64 = 0.0 # [0, 1]
essentialProteinProductionRate::Float64 = 0.0
essentialProteinDegradationRate::Float64 = 0.0 # 1/h

# essential metabolite kinetics
essentialMetaboliteProductionRate::Float64 = 2.0
essentialMetaboliteDegradationRate::Float64 = 0.5
essentialMetaboliteKm::Float64 = 50.0
essentialMetaboliteThreshold::Float64 = 0.0

# simulation settings
agentTimeStep::Float64 = 5/60 # h
duration = 10.0
showProgress::Bool = true

# Init parameter struct
parameters::BioreactorParameters = BioreactorParameters(
    startingVolume, startingCellCount, startingEssentialProteinConcentration,
    muMax, carryingCapacity,
    divisionSymmetry, essentialProteinProductionRate, essentialProteinDegradationRate,
    essentialMetaboliteProductionRate, essentialMetaboliteDegradationRate, essentialMetaboliteKm, essentialMetaboliteThreshold,
    agentTimeStep, duration, showProgress
)

# Init Bioreactor struct
b = Bioreactor(parameters)
# Run bioreactor
simulateBioreactor(b)

# Plotting
include("./Plotting/plotting.jl")
#display(plotBiomass(b))
display(plotBioreactor(b))
#display(essentialProteinHistogram(b, 3.0))

#=
timepoints = collect(0:5/60:b.parameters.duration)
biomass = [sum(b.solution(t)[getCellIdx():3:end]) for t in timepoints]
dBdt = calculateDerivative(timepoints, biomass)
dbdt = calculateDerivative(timepoints, dBdt)
plot(timepoints, biomass)
plot(timepoints, dBdt)
plot(timepoints, dbdt)
=#

k = 10
timepoints = collect(0:5/60:b.parameters.duration)
biomass = [sum(b.solution(t)[getCellIdx():3:end]) for t in timepoints]
timepoints = [mean(timepoints[(i-k):(i+k)]) for i in (1+k):(length(timepoints)-k)]
biomass = [mean(biomass[(i-k):(i+k)]) for i in (1+k):(length(biomass)-k)]
biomass ./= maximum(biomass)
dBdt = calculateDerivative(timepoints, biomass)
dbdt = calculateDerivative(timepoints, dBdt)
display(plot(timepoints, biomass))
display(plot(timepoints, dBdt))
display(plot(timepoints, dbdt))