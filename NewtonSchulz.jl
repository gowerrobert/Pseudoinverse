## NEWTON SCHULZ
function boot_NewtonSchulz(prob::Prob,options::MyOptions)
    if(issparse(prob.A))
        flopsperiter = nnz(prob.A)*options.m + options.m^2*options.n+2*options.n*options.m;#nnz(A)*m +m^2*n +2*n*m;
    else
        flopsperiter =  2*options.m^2*options.n+2*options.n*options.m;#nnz(A)*m +m^2*n +2*n*m;
    end
    M0 = (prob.A')/(2*trace(prob.A*prob.A'));
    method= Method( flopsperiter, "NewtonSchulz",M0, iter_NewtonSchulz);
    options.skip_error_calculation =1; # No skipping error calculations since iterations are already so slow
    return method
end
function iter_NewtonSchulz(A,M,options,iteration)
    M =2*M - M*A*M;
end