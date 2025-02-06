using Plots
include("./bioreactor.jl")

# initial conditions
startingCellCount::Int64 = 100 # cells
startingVolume::Float64 = 400*10^-11 # L
startingEssentialProteinConcentration::Float64 = 800 # molecules

# growth kinetics
muMax::Float64 = 2 # divisions/h
carryingCapacity::Float64 = 1e6

# essentialProtein kinetics
divisionSymmetry::Float64 = 1.0 # [0, 1]
essentialProteinDegradationRate::Float64 = 0.05 # 1/h

# essential metabolite kinetics
essentialMetaboliteProductionRate::Float64 = 4.0
essentialMetaboliteDegradationRate::Float64 = 2.0
essentialMetaboliteKm::Float64 = 50

# simulation settings
agentTimeStep::Float64 = 5/60 # h
duration = 3.0
showProgress::Bool = true

parameters::BioreactorParameters = BioreactorParameters(
    startingVolume, startingCellCount, startingEssentialProteinConcentration,
    muMax, carryingCapacity,
    divisionSymmetry, essentialProteinDegradationRate,
    essentialMetaboliteProductionRate, essentialMetaboliteDegradationRate, essentialMetaboliteKm,
    agentTimeStep, duration, showProgress
)


b = Bioreactor(parameters)
simulateBioreactor(b, duration)

include("./Plotting/plotting.jl")
#display(plotBiomass(b))
display(plotBioreactor(b))
display( essentialProteinHistogram(b, duration))
savefig(plotBioreactor(b), "./Figures/temp.png")
