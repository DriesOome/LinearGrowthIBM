# single plots
function plotRonald(b1::Bioreactor, b2::Bioreactor)
    fontsize=14
    p = plot(layout=(1,2), framestyle=:box)
    ylim_cellDensity = sum(b1.solution(b1.parameters.duration)[getCellIdx():3:end])

    # mappings
    time = collect(0:b1.parameters.agentTimeStep:b1.parameters.duration)
    p1 = plot(time, [sum(b1.solution(t)[getCellIdx():3:end]) for t in time], c=:blue, xlabel="time (h)", ylabel="density/L", title="exponential growth",
        xticks=[i for i in 0:3],
        ylims=[0, ylim_cellDensity], legend=:none,left_margin = 5mm,  bottom_margin = 5mm)
    yaxis2 = twinx()
    p1 = plot!(yaxis2, time, [100*getFractionGrowingCells(b1.solution(t), b1.parameters.essentialMetaboliteThreshold) for t in time],
        ylabel="", ylims=[0,105], c=:red, yformatter=_->"", legend=:none)
    
    time = collect(0:b2.parameters.agentTimeStep:b2.parameters.duration)
    p2 = plot(time, [sum(b2.solution(t)[getCellIdx():3:end]) for t in time], c=:blue, xlabel="time (h)", ylabel="", ylims=[0, ylim_cellDensity], title="linear growth",
        xticks=[i for i in 0:2:15],
        yformatter=_->"", legend=:none)
    yaxis2 = twinx()
    p2 = plot!(yaxis2, time, [100*getFractionGrowingCells(b2.solution(t), b2.parameters.essentialMetaboliteThreshold) for t in time],
        ylabel="% growing", ylims=[0,105], c=:red, legend=:none, right_margin = 5mm)

    p = plot!(p1, p2, size=(800,400), dpi=1000)
    #p1 = plot!(p1, xguidefontsize=fontsize, yguidefontsize=fontsize, legendfontpointsize=fontsize, legend=:none)
    #p1 = plot!(p1, xtick=[i for i in 0:1:100])
    return p
end
