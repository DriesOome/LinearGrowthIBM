using Random
using ProgressBars
using DifferentialEquations
using Distributions
include("./bioreactorParameters.jl")

mutable struct Bioreactor
    # time::Float64
    ## reactor variables
    # substrateConcentration::Float64
    # volume::Float64
    # productConcentration::Float64
    ## cell variables
    # volume::Float64
    # essentialProtein::Float64

    parameters::BioreactorParameters
    solution
    progressBar::ProgressBar

    function Bioreactor(parameters::BioreactorParameters)
        bioreactor::Bioreactor = new(parameters, undef,
        ProgressBar(total=Int64(parameters.duration/parameters.agentTimeStep)))
        return bioreactor
    end
end

# init functions
function initializeBacterialPopulation(bioreactor::Bioreactor)
    u_bacteria::Vector{Float64} = []
    for i in 1:bioreactor.parameters.startingCellCount
        push!(u_bacteria, rand()+1) # start with random volume
        push!(u_bacteria, bioreactor.parameters.startingEssentialProteinConcentration) # starting essential protein count
        push!(u_bacteria, bioreactor.parameters.essentialMetaboliteProductionRate*bioreactor.parameters.startingEssentialProteinConcentration/(bioreactor.parameters.essentialMetaboliteDegradationRate)) # starting essential metabolite count
    end
    return u_bacteria
end

# Function defining the initial state of the bioreactor
function initializeBioreactorState(bioreactor::Bioreactor)
    u::Vector{Float64} = []
    # global state variables
    push!(u, bioreactor.parameters.startingVolume) # volume
    push!(u, 0) # product
    # cell state variables
    append!(u, initializeBacterialPopulation(bioreactor))
end

function simulateBioreactor(bioreactor::Bioreactor, duration::Float64)
    # inital state
    u_init = initializeBioreactorState(bioreactor)
    # setup agent actions https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/
    stops = collect(0:bioreactor.parameters.agentTimeStep:duration)
    agentActionCondition(u, t, integrator) = t in stops
    agentActionAffect!(integrator) = agentActions(integrator, bioreactor)
    cb = DiscreteCallback(agentActionCondition, agentActionAffect!;
        save_positions = (true, true))
    bioreactor.solution = solve(ODEProblem(bioreactorODEFunction, u_init, (0, duration), bioreactor),
        callback = cb, tstops = stops
    )
end

function agentActions(integrator, bioreactor::Bioreactor)
    divideCells(integrator, bioreactor)
    # update the progress bar
    if bioreactor.parameters.showProgress == true; update(bioreactor.progressBar) end
    nothing
end

function divideCells(integrator, bioreactor)
    # all cells with a volume larger than 1 are dividing
    dividingCells = findall(integrator.u[getCellIdx():3:end] .> 2)
    for parentId in dividingCells
        resize!(integrator, length(integrator.u)+3)
        childId = totalCells(integrator.u)
        # reset volumes
        setCellVolume!(integrator.u, childId,  getCellVolume(integrator.u, parentId)-rand(Normal(1,0.2)))
        setCellVolume!(integrator.u, parentId, getCellVolume(integrator.u, parentId)-getCellVolume(integrator.u, childId))
        # divide essential protein         
        setCellEssentialProtein!(integrator.u, childId, rand(Binomial(floor(getCellEssentialProtein(integrator.u, parentId)), bioreactor.parameters.divisionSymmetry)))
        setCellEssentialProtein!(integrator.u, parentId, getCellEssentialProtein(integrator.u, parentId)-getCellEssentialProtein(integrator.u, childId))
        # set metabolite
        setCellEssentialMetabolite!(integrator.u, childId, getCellEssentialMetabolite(integrator.u, parentId)*0.5)
        setCellEssentialMetabolite!(integrator.u, parentId, getCellEssentialMetabolite(integrator.u, parentId)*0.5)
    end
    return nothing
end

function bioreactorODEFunction(du, u, bioreactor::Bioreactor, t)
    # simulate bioreactor variables
    du[getVolumeIdx()] = calculateChangeVolume(bioreactor)
    du[getProductIdx()] = calculateChangeProductConcentration(bioreactor)
    # simulate cell variables
    essentialProteinCounts::Vector{Float64} = u[getCellIdx()+1:3:end]
    essentialMetaboliteCounts::Vector{Float64} = u[getCellIdx()+2:3:end]
    essentialMetaboliteConcentration::Vector{Float64} = essentialMetaboliteCounts./u[getCellIdx():3:end]
    growthModifications::Vector{Float64} = broadcast(max, 0.0, (essentialMetaboliteConcentration)./((essentialMetaboliteConcentration) .+ bioreactor.parameters.essentialMetaboliteKm))
    thresholdBits::Vector{Float64} = essentialMetaboliteConcentration .>= bioreactor.parameters.essentialMetaboliteThreshold
    du[getCellIdx():3:end] = bioreactor.parameters.muMax.*growthModifications.*thresholdBits
    du[getCellIdx()+1:3:end] = -bioreactor.parameters.essentialProteinDegradationRate.*essentialProteinCounts
    du[getCellIdx()+2:3:end] = bioreactor.parameters.essentialMetaboliteProductionRate.*essentialProteinCounts-bioreactor.parameters.essentialMetaboliteDegradationRate.*essentialMetaboliteCounts
    nothing
end 


# ODE function
function calculateChangeVolume(bioreactor::Bioreactor)
    dV::Float64 = 0.0
    return dV
end

function calculateChangeProductConcentration(bioreactor::Bioreactor)
    dP::Float64 = 0.0
    return dP
end
# helper functions

# getters and setters
function totalCells(u)
    return Int64((length(u)-getProductIdx())/3)
end

function getCellIdx(cellId)
    return Int64(3*(cellId-1) + getCellIdx())
end

function getCellVolume(u, cellId)
    return u[getCellIdx(cellId)]
end

function setCellVolume!(u, cellId, volume)
    u[getCellIdx(cellId)] = volume
end

function getCellEssentialProtein(u, cellId)
    return u[getCellIdx(cellId)+1]
end

function setCellEssentialProtein!(u, cellId, essentialProtein)
    u[getCellIdx(cellId)+1] = essentialProtein
end

function getCellEssentialMetabolite(u, cellId)
    return u[getCellIdx(cellId)+2]
end

function setCellEssentialMetabolite!(u, cellId, essentialMetabolite)
    u[getCellIdx(cellId)+2] = essentialMetabolite
end

function getVolumeIdx()
    return 1
end

function getProductIdx()
    return getVolumeIdx()+1
end

function getCellIdx()
    return getProductIdx()+1
end

function getFractionGrowingCells(u, growthThreshold::Float64)
    fraction = 0
    for cellId in 1:totalCells(u)
        if getCellEssentialMetabolite(u, cellId)/getCellVolume(u, cellId) >= growthThreshold
            fraction += 1
        end
    end
    fraction = fraction/totalCells(u)
    return fraction
end