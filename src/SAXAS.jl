## THE SAXAS method for symmetric matrices
function boot_SAXAS(prob,options)
#    options = set_quNac_standard_options(A,options);
    options.sketchsize = ceil(sqrt(options.n));
    flopsperiter = 7*options.n^2*options.sketchsize +2*options.n*options.sketchsize^2 +options.sketchsize^3+3*options.n^2+options.n;
    name = "SAXAS-$(options.sketch)";
    if(options.M0type == "AAproj")
         M0 = (trace(prob.A)*(prob.A)*(prob.A))/vecnorm(prob.A)^4;  
         name = "$(name)-$(options.M0type)"
    elseif (options.M0type == "Aproj")
         M0 = (min(options.n,options.m)*prob.A')/vecnorm(prob.A)^2;
         name = "$(name)-$(options.M0type)"
    elseif(options.M0type == "NS")
         M0 =  (prob.A')/(2*trace(prob.A*prob.A'));
         name = "$(name)-$(options.M0type)"
    elseif (options.M0type == "eye")
         M0 = eye(prob.A');
         name = "$(name)-$(options.M0type)"        
    elseif (options.M0type == "ATAnorm")
         M0 = (prob.A)*(prob.A)/vecnorm(prob.A)^2;  
        # name = "$(name)-$(options.M0type)"  
    else
         println("Unknown type of M0!!")
    end
    method= Method(flopsperiter, name,M0, iter_SAXAS) ;
    return method;
end

# The SAXAS iteration!
function iter_SAXAS(A,M, options,iteration)
    # Select which sketch to use;
    local s::Array{Int64}=sample(1:options.n,options.sketchsize,replace=false);    
    if(options.sketch =="uni" )
        AS = A[:,s]; 
        elseif(options.sketch =="ada" )
        AS = A*M[:,s];  R = M[:,s];
        #elseif(options.sketch == "MT" )   # Same thing in symmetric case!
         #       AS = A*(M[s,:]');  R = M[s,:]';
else
        println("None existant sampling in iter_SAXAS! Exiting");
    return;
end
    if(options.sketch =="uni")
        SAS = A[s,s];
    else
        SAS = (AS')*R;
    end
SAAS = AS'*AS;
    if(issparse(A))
    invSAAS =  pinv(full(SAAS));
 else
     invSAAS =  pinv(SAAS);
 end
ASinvSAAS = AS*invSAAS;
SAMAS = AS'*(M*AS);
# M = M + RinvA*(R') - RinvA*RAMAR*(RinvA'); % MAR*(RinvA')+RinvA*((MAR')*Aproj);
## Testing - Comparing
#M = M + RinvA*(R')*(A -A*M*A )*(RinvA*(R')); % MAR*(RinvA')+RinvA*((MAR')*Aproj);
#M = M + RinvA*(AR' -AR'*M*A )*(RinvA*(R'));
M = M + ASinvSAAS*(SAS -SAMAS )*(ASinvSAAS');
end
