include("bioreactor.jl")
using DataFrames
using CSV

function extractDataFrame(bioreactor::Bioreactor, timePoints::Vector{Float64})
    data::DataFrame = DataFrame(time=[], cellId=[], volume=[], essentialProtein=[], essentialMetabolite=[])
    for timePoint in timePoints
        u = bioreactor.solution(timePoint)
        for cellId in 1:totalCells(u)
            push!(data, [timePoint, cellId, getCellVolume(u, cellId), getCellEssentialProtein(u, cellId), getCellEssentialMetabolite(u, cellId)])
        end
    end
    return data
end

function saveDataFrame(data::DataFrame, filename::String)
    CSV.write(filename, data)
end

function saveParameters(p::BioreactorParameters, filename)
    CSV.write(filename, DataFrame(
        startingVolume = p.startingVolume,
        startingCellCount = p.startingCellCount,
        startingEssentialProteinConcentration = p.startingEssentialProteinConcentration,
        
        muMax = p.muMax,
        carryingCapacity = p.carryingCapacity, 

        divisionSymmetry = p.divisionSymmetry,
        essentialProteinDegradationRate = p.essentialProteinDegradationRate,

        essentialMetaboliteProductionRate = p.essentialMetaboliteProductionRate,
        essentialMetaboliteDegradationRate = p.essentialMetaboliteDegradationRate,
        essentialMetaboliteKm = p.essentialMetaboliteKm,
        essentialMetaboliteThreshold = p.essentialMetaboliteThreshold,

        agentTimeStep = p.agentTimeStep,
        duration = p.duration,
    ))
end

function saveBioreactorRun(bioreactor::Bioreactor, timePoints::Vector{Float64}, filename::String)
    data = extractDataFrame(bioreactor, timePoints)
    saveDataFrame(data, filename*"_data.csv")
    saveParameters(bioreactor.parameters, filename*"_params.csv")
end