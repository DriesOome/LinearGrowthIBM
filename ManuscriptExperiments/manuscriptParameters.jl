include("../bioreactorParameters.jl")
# initial conditions
startingCellCount::Int64 = 100 # cells
startingVolume::Float64 = 400*10^-11 # L
startingEssentialProteinConcentration::Float64 = 200 # molecules

# growth kinetics
muMax::Float64 = 2 # divisions/h
carryingCapacity::Float64 = 1e6

# essentialProtein kinetics
divisionSymmetry::Float64 = 0.5 # [0, 1]
essentialProteinProductionRate::Float64 = 0.0
essentialProteinDegradationRate::Float64 = 0.05 # 1/h

# essential metabolite kinetics
essentialMetaboliteProductionRate::Float64 = 2.0
essentialMetaboliteDegradationRate::Float64 = 0.5
essentialMetaboliteKm::Float64 = 50
essentialMetaboliteThreshold::Float64 = essentialMetaboliteKm/4

# simulation settings
agentTimeStep::Float64 = 5/60 # h
duration = 15.0
showProgress::Bool = true

manuscriptParameters::BioreactorParameters = BioreactorParameters(
    startingVolume, startingCellCount, startingEssentialProteinConcentration,
    muMax, carryingCapacity,
    divisionSymmetry, essentialProteinProductionRate, essentialProteinDegradationRate,
    essentialMetaboliteProductionRate, essentialMetaboliteDegradationRate, essentialMetaboliteKm, essentialMetaboliteThreshold,
    agentTimeStep, duration, showProgress
)