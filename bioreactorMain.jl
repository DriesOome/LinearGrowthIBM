using Plots
include("./bioreactor.jl")

# initial conditions
startingCellCount::Int64 = 1000 # cells
startingVolume::Float64 = 400*10^-11 # L
startingEssentialProteinConcentration::Float64 = 200 # molecules

# growth kinetics
muMax::Float64 = 2.0#log(2) # divisions/h
carryingCapacity::Float64 = 1e8

# essentialProtein kinetics
divisionSymmetry::Float64 = 0.5 # [0, 1]
essentialProteinProductionRate::Float64 = 0.0
essentialProteinDegradationRate::Float64 = 0.0 # 1/h

# essential metabolite kinetics
essentialMetaboliteProductionRate::Float64 = 2.0
essentialMetaboliteDegradationRate::Float64 = 0.5
essentialMetaboliteKm::Float64 = 50.0
essentialMetaboliteThreshold::Float64 = 5

# simulation settings
agentTimeStep::Float64 = 5/60 # h
duration = 30.0
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
display(plotBioreactor(b))


