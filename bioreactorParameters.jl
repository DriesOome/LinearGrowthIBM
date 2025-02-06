mutable struct BioreactorParameters
    # initial conditions
    startingVolume::Float64
    startingCellCount::Int64
    startingEssentialProteinConcentration::Float64

    # growth kinetics
    muMax::Float64 # divisions/h
    carryingCapacity::Float64 

    # essentialProtein kinetics
    divisionSymmetry::Float64
    essentialProteinDegradationRate::Float64

    # essential metabolite kinetics
    essentialMetaboliteProductionRate::Float64
    essentialMetaboliteDegradationRate::Float64
    essentialMetaboliteKm::Float64
    essentialMetaboliteThreshold::Float64

    # simulation settings
    agentTimeStep::Float64 # h
    duration::Float64 # h
    showProgress::Bool
end

Base.copy(p::BioreactorParameters) = 
    BioreactorParameters(
        # initial conditions
        p.startingVolume, p.startingCellCount, p.startingEssentialProteinConcentration, 
        # growth kinetics
        p.muMax, p.carryingCapacity,
        # essentialProtein kinetics
        p.divisionSymmetry, p.essentialProteinDegradationRate,
        # essential metabolite kinetics
        p.essentialMetaboliteProductionRate, p.essentialMetaboliteDegradationRate, p.essentialMetaboliteKm, essentialMetaboliteThreshold,
        # simulation settings
        p.agentTimeStep, p.duration, p.showProgress
    )