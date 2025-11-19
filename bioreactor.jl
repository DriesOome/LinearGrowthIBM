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
function initializeBacterialPopulation(b::Bioreactor)
    u_bacteria::Vector{Float64} = []
    essentialMetabolite::Float64 = solve(NonlinearProblem((u, p) -> p[1]*p[2] .- p[3].*u - p[4].*(u./(u .+ p[5])), [0.0], [
        b.parameters.essentialMetaboliteProductionRate, b.parameters.startingEssentialProteinConcentration, b.parameters.essentialMetaboliteDegradationRate, b.parameters.muMax, b.parameters.essentialMetaboliteKm
    ]))[1]
    for i in 1:b.parameters.startingCellCount
        cellVolume::Float64 = rand(Uniform(0,1))+1
        push!(u_bacteria, cellVolume) # start with random volume
        push!(u_bacteria, b.parameters.startingEssentialProteinConcentration) # starting essential protein count
        push!(u_bacteria, essentialMetabolite*cellVolume) # starting essential metabolite count
    end
    return u_bacteria
end

# Function defining the initial state of the bioreactor
function initializeBioreactorState(bioreactor::Bioreactor)
    u::Vector{Float64} = []
    # global state variables
    push!(u, bioreactor.parameters.startingVolume) # volume
    push!(u, 0) # dead cell count
    push!(u, 0) # dead cell volume
    # cell state variables
    append!(u, initializeBacterialPopulation(bioreactor))
end

function simulateBioreactor(bioreactor::Bioreactor)
    # inital state
    u_init = initializeBioreactorState(bioreactor)
    # setup agent actions https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/
    stops = collect(0:bioreactor.parameters.agentTimeStep:bioreactor.parameters.duration)

    agentActionCondition(u, t, integrator) = t in stops
    agentActionAffect!(integrator) = agentActions(integrator, bioreactor)
    agentCallback = DiscreteCallback(agentActionCondition, agentActionAffect!;
        save_positions = (true, true))

    stoppingCondition(u, t, integrator) = totalCells(u) > 1e5
    stoppingAffect!(integrator) = terminateBioreactor(integrator, bioreactor)
    stoppingCallback = DiscreteCallback(stoppingCondition, stoppingAffect!)

    prob = ODEProblem(bioreactorODEFunction, u_init, (0, bioreactor.parameters.duration), bioreactor)
    integrator = init(prob, Tsit5(), tstops=stops)

    ts = Float64[]
    us = Vector{Vector{Float64}}()
    push!(ts, integrator.t)
    push!(us, copy(integrator.u))
    while integrator.t < bioreactor.parameters.duration
        step!(integrator)
        divideCells(integrator, bioreactor)
        removeDeadCells(integrator, bioreactor)
        u_modified!(integrator, true)
        push!(ts, integrator.t)
        push!(us, copy(integrator.u))
    end
    bioreactor.solution = DiffEqBase.build_solution(prob, integrator.alg, ts, us)
end

function terminateBioreactor(integrator, bioreactor::Bioreactor)
    bioreactor.parameters.duration = integrator.t
    terminate!(integrator)
end

function agentActions(integrator, bioreactor::Bioreactor)
    divideCells(integrator, bioreactor)
    # update the progress bar
    if bioreactor.parameters.showProgress == true; update(bioreactor.progressBar) end
    nothing
end


function removeDeadCells(integrator, bioreactor::Bioreactor)
    u = integrator.u
    essentialMetaboliteCounts::Vector{Float64} = u[getCellIdx()+2:3:end]
    essentialMetaboliteConcentration::Vector{Float64} = essentialMetaboliteCounts./u[getCellIdx():3:end]
    growthModifications::Vector{Float64} = broadcast(max, 0.0, (essentialMetaboliteConcentration)./(essentialMetaboliteConcentration .+ bioreactor.parameters.essentialMetaboliteKm))
    nonGrowingCellsIds = findall(growthModifications .< 0.01)
    integrator.u[getDeadCellCountIdx()] += length(nonGrowingCellsIds)
    integrator.u[getDeadCellVolumeIdx()] += sum([getCellVolume(integrator.u, cellId) for cellId in nonGrowingCellsIds])
    nonGrowingCellsIds = nonGrowingCellsIds .- 1
    nonGrowingCellsIds = nonGrowingCellsIds.*3 .+ getCellIdx()
    toDelete = copy(nonGrowingCellsIds)
    append!(toDelete, nonGrowingCellsIds .+ 1)
    append!(toDelete, nonGrowingCellsIds .+ 2)
    toDelete = sort(toDelete)
    deleteat!(integrator, toDelete)
end

function divideCells(integrator, bioreactor)
    # all cells with a volume larger than 1 are dividing
    dividingCells = findall(integrator.u[getCellIdx():3:end] .> 2)
    for parentId in dividingCells
        resize!(integrator, length(integrator.u)+3)
        childId = totalGrowingCells(integrator.u)
        # reset volumes
        setCellVolume!(integrator.u, childId,  getCellVolume(integrator.u, parentId)-rand(Normal(1,0.05)))
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
    du[getDeadCellCountIdx()] = 0.0
    du[getDeadCellVolumeIdx()] = 0.0
    # simulate cell variables
    essentialProteinCounts::Vector{Float64} = u[getCellIdx()+1:3:end]
    essentialMetaboliteCounts::Vector{Float64} = u[getCellIdx()+2:3:end]
    essentialMetaboliteConcentration::Vector{Float64} = essentialMetaboliteCounts./u[getCellIdx():3:end]
    growthModifications::Vector{Float64} = broadcast(max, 0.0, (essentialMetaboliteConcentration)./(essentialMetaboliteConcentration .+ bioreactor.parameters.essentialMetaboliteKm))
    thresholdBits::Vector{Float64} = essentialMetaboliteConcentration .>= bioreactor.parameters.essentialMetaboliteThreshold
    du[getCellIdx():3:end] = (bioreactor.parameters.muMax.*growthModifications.*thresholdBits./log(2))
    du[getCellIdx()+1:3:end] = bioreactor.parameters.essentialProteinProductionRate .- bioreactor.parameters.essentialProteinDegradationRate.*essentialProteinCounts
    du[getCellIdx()+2:3:end] = bioreactor.parameters.essentialMetaboliteProductionRate.*essentialProteinCounts .- bioreactor.parameters.essentialMetaboliteDegradationRate.*essentialMetaboliteCounts
    nothing
end 


# ODE function
function calculateChangeVolume(bioreactor::Bioreactor)
    dV::Float64 = 0.0
    return dV
end


# helper functions

# getters and setters
function totalGrowingCells(u)
    return Int64((length(u)-getDeadCellVolumeIdx())/3)
end

function totalCells(u)
    return Int64((length(u)-getDeadCellVolumeIdx())/3) + u[getDeadCellCountIdx()]
end

function getTotalCellDensity(u)
    return sum(u[getCellIdx():3:end]) + u[getDeadCellVolumeIdx()]
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

function getDeadCellCountIdx()
    return getVolumeIdx()+1
end

function getDeadCellVolumeIdx()
    return getDeadCellCountIdx()+1
end

function getCellIdx()
    return getDeadCellVolumeIdx()+1
end

function getFractionGrowingCells(u, growthThreshold::Float64)
    return totalGrowingCells(u)/totalCells(u)
end