using Plots
using DifferentialEquations
using Distributions
using Base.Threads
include("./bioreactor.jl")

# initial conditions
startingCellCount::Int64 = 100 # cells
startingVolume::Float64 = 400*10^-11 # L
startingEssentialProteinConcentration::Float64 = 20 # molecules

# growth kinetics
muMax::Float64 = 2 # divisions/h
carryingCapacity::Float64 = 1e6

# essentialProtein kinetics
divisionSymmetry::Float64 = 0.0 # [0, 1]
essentialProteinDegradationRate::Float64 = 0.0 # 1/h

# essential metabolite kinetics
essentialMetaboliteProductionRate::Float64 = 10.0
essentialMetaboliteDegradationRate::Float64 = 1.0
essentialMetaboliteKm::Float64 = 10

# simulation settings
agentTimeStep::Float64 = 5/60 # h
duration = 20.0
showProgress::Bool = false

parameters::BioreactorParameters = BioreactorParameters(
    startingVolume, startingCellCount, startingEssentialProteinConcentration,
    muMax, carryingCapacity,
    divisionSymmetry, essentialProteinDegradationRate,
    essentialMetaboliteProductionRate, essentialMetaboliteDegradationRate, essentialMetaboliteKm,
    agentTimeStep, duration, showProgress
)


variableParameter::Vector{Float64} = [0.0,0.5]

bioreactors::Vector{Bioreactor} = []
for i in eachindex(variableParameter)
    p = copy(parameters)
    p.divisionSymmetry = variableParameter[i]
    push!(bioreactors, Bioreactor(p))
end
progress = ProgressBar(total=length(variableParameter))
Threads.@threads for i in eachindex(variableParameter)
    simulateBioreactor(bioreactors[i], duration)
    update(progress)
end


include("./Plotting/plotting.jl")
display(plotBioreactors(bioreactors))

