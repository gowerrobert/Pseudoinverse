# Front end for plotting the execution in time and in flops of the outputs recorded in OUTPUTS.
function plot_outputs_Plots(OUTPUTS,probname)
    output = OUTPUTS[1];
    plot(output.times',output.residuals', 
    xlabel = "time", 
    ylabel = "residual", 
    yscale = :log10,
    label  = output.name,
    linestyle=:auto,   tickfont=font(14), guidefont=font(18),
    marker =:auto, 
    grid = false)
    for i =2:length(OUTPUTS)
        output = OUTPUTS[i]; 
        plot!(output.times',output.residuals', yscale = :log10, label  = output.name, linestyle=:auto, marker =:auto, grid = false)
    end
    println(probname)
    probname= replace(probname, r"[\/]", "-");
    savefig("../figures/$(probname).pdf");
    # Now in residual X flops
    output = OUTPUTS[1];
    lt = length(output.times);
    reflops = output.flopsperiter;
    plot((output.flopsperiter/reflops)*(1:lt)*(output.iterations/lt),output.residuals', 
    xlabel = "flops", 
    ylabel = "residual", 
    yscale = :log10,
    label  = output.name,
    linestyle=:auto,  tickfont=font(14), guidefont=font(18),
    marker =:auto, 
    grid = false)
    for i =2:length(OUTPUTS)
        output = OUTPUTS[i];
        lt = length(output.times);
        plot!((output.flopsperiter/reflops)*(1:lt)*(output.iterations/lt),output.residuals',
        yscale = :log10, label  = output.name, linestyle=:auto, marker =:auto, grid = false)
    end
#    probname= replace(probname, r"[\/]", "-");
    savefig("../figures/$(probname)-flops.pdf");
    
end


