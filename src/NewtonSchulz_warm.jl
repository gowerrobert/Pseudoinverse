## NEWTON SCHULZ with Randomized SATAX warm starting
function boot_NewtonSchulz_warm(prob::Prob,options::MyOptions)
 
    method =boot_SATAX(prob,options);
    method.stepmethod = iter_NewtonSchulz_warm;
    method.name = "NS-$(method.name)" 
    return method
end
function iter_NewtonSchulz_warm(A,M,options,iteration)
    
    if(options.sketch == "uni")
        datapass = options.sketchsize/options.n
    elseif(options.sketch == "ada")
        datapass = options.sketchsize/options.m
    end
    println("datapass percentage: $(datapass*iteration)")
    if(datapass*iteration >2.0)
       #println("NS -step")
       M =2*M - M*A*M;
    else
       #println("SATAX -step")
       M = iter_SATAX(A,M, options,iteration);
    end
end