# single plots
function plotManuscriptFigure(b1::Bioreactor, b2::Bioreactor)
    fontsize=14
    p = plot(framestyle=:box)
    ylim_cellDensity = sum(b1.solution(b1.parameters.duration)[getCellIdx():3:end])

    # mappings
    time = collect(0:b1.parameters.agentTimeStep:b1.parameters.duration)
    p1 = plot(time, [sum(b1.solution(t)[getCellIdx():3:end]) for t in time], c=:blue, 
        xlabel="Time (h)", ylabel="density/L", title="exponential growth", legend=:none,
        xguidefontsize=14, yguidefontsize=14, titlefontsize=16,
        xticks=[i for i in 0:3], ylims=[0, ylim_cellDensity], 
        left_margin = 5mm,  bottom_margin = 5mm)
    yaxis2 = twinx()
    p1 = plot!(yaxis2, time, [100*getFractionGrowingCells(b1.solution(t), b1.parameters.essentialMetaboliteThreshold) for t in time],
        ylabel="% growing", yguidefontsize=14,  ylims=[0,105], c=:red, legend=:none)

    time = collect(0:b2.parameters.agentTimeStep:b2.parameters.duration)
    p2 = plot(time, [sum(b2.solution(t)[getCellIdx():3:end]) for t in time], c=:blue, 
        xlabel="Time (h)", ylabel="density/L", title="linear growth",
        xguidefontsize=14, yguidefontsize=14, titlefontsize=16,
        xticks=[i for i in 0:2:15], legend=:none, ylims=[0, ylim_cellDensity])
    yaxis2 = twinx()
    p2 = plot!(yaxis2, time, [100*getFractionGrowingCells(b2.solution(t), b2.parameters.essentialMetaboliteThreshold) for t in time],
        ylabel="% growing", yguidefontsize=14, ylims=[0,105], c=:red, legend=:none, right_margin = 5mm)

    p = plot!(p1, p2, size=(500,800), dpi=1000, layout=(2,1))
    #p1 = plot!(p1, xguidefontsize=fontsize, yguidefontsize=fontsize, legendfontpointsize=fontsize, legend=:none)
    #p1 = plot!(p1, xtick=[i for i in 0:1:100])
    return p
end
