## THE SATAX method
function boot_SATAX(prob,options)
    if(options.sketch == "uni")
        options.sketchsize = ceil(sqrt(options.n));  # sample over columns of A
        flopsperiter = 4*options.n*options.m *options.sketchsize+
        2*options.n*options.sketchsize^2+options.sketchsize^3+options.n*options.m;
    elseif(options.sketch == "ada")
        options.sketchsize = ceil(sqrt(options.m));  # sample over m columns of M
            flopsperiter = 5*options.n*options.m *options.sketchsize+
        2*options.n*options.sketchsize^2+options.sketchsize^3+options.n*options.m;
        
    elseif(options.sketch == "MMT" )
        options.sketchsize = ceil(sqrt(options.n));  # sample over rows of M
            flopsperiter = 6*options.n*options.m *options.sketchsize+
        2*options.n*options.sketchsize^2+options.sketchsize^3+options.n*options.m; 
    else
        println("None existant sampling in iter_Rightside! Exiting");
        return;
    end
 # 4 mXtauXn + mXn + 2tau^2Xn+tauXm + tau^3
    name = "SATAX-$(options.sketch)";
    if(options.M0type == "ATnorm")
         M0 =  prob.A'/vecnorm(prob.A);   
         name = "$(name)-$(options.M0type)";
    elseif(options.M0type == "NS")
         M0 =  (prob.A')/(2*trace(prob.A*prob.A'));
         name = "$(name)-$(options.M0type)"
    #elseif(options.M0type == "ATAproj") #NOTE POSSIBLE! Dimension mismatch
    #     M0 =  ((prob.A')*prob.A)/vecnorm(prob.A)^2;
    #     name = "$(name)-$(options.M0type)"
    elseif(options.M0type == "ATproj")    
       M0 = (min(options.n,options.m)*prob.A')/vecnorm(prob.A)^2; 
      # name = "$(name)-$(options.M0type)"
    else
        println("unknown type of M0!!")
    end

    method= Method(flopsperiter, name,M0, iter_SATAX) ;
    return method;
end

# X_{k+1} =X_k-A^\top A S (S^\top A^\top A A^\top AS)^{\dagger}S^\top A^\top(AX_k-I).
function iter_SATAX(A,M, options,iteration)
    # Select which sketch to use;
    local s::Array{Int64}=[1];
    if(options.sketch =="uni" )
        s = sample(1:options.n,options.sketchsize,replace=false);  # sample over columns
        AS = A[:,s]; 
    elseif(options.sketch =="ada" )
        s= sample(1:options.m,options.sketchsize,replace=false);  # sample over m columns of X
        AS = A*M[:,s]; 
elseif(options.sketch == "XXT" )
        s= sample(1:options.n,options.sketchsize,replace=false); # sample over columns
        AS = A*(M*M[s,:]');     
else
        println("None existant sampling in iter_Rightside! Exiting");
    return;
end
R = (AS'*A)*M-AS';   #S^\top A^\top(AX_k-I)  # 2mXtauXn +tauXm
ATSTA= AS'*A; # A'*AS);% Think about avoiding transposing.  # mXtauXn
STATAATAS = (ATSTA*ATSTA');  # n X tau^2
invSTATAATAS =  pinv(full(STATAATAS)); # tatu^3
M = M - (ATSTA')*invSTATAATAS*R; # tau^2 * n +  nXmXtau + mXn
# Testing - Comparing
# M = M + RinvA*(R')*(A -A*M*A )*(RinvA*(R')); % MAR*(RinvA')+RinvA*((MAR')*Aproj);
    #total flops:  4 mXtauXn + mXn + 2tau^2Xn+tauXm + tau^3
end
