# A wrapper function for testing and timing iterative methods for
# calculating the pseudoinverse of a matrix - 2016 - Robert M. Gower
# PseudoInvRand Copyright (C) 2016, Robert Gower
function  pseudoinvert_matrix(prob::Prob, method_name, options::MyOptions )

    method = boot_method(method_name,prob,options); 
    if(method=="METHOD DOES NOT EXIST")
        println("FAIL: unknown method name:")
        return
    end
    local M= method.M0;
    println(method.name);
#     times[1] = toc;
    times= [0];
    if(options.exacterror)
        initial_error= vecnorm(prob.Apseudo-M);
    end
    initial_residual= vecnorm(prob.A -prob.A*M*prob.A);
    errors = [1];
    residuals  = [1];
    local tickcounter =1;
    local timeaccum=0;
    iterations =0;
    fail = "failed";
# # Print heard
    if(options.printiters)
         println("-------------------")
         println("It   | Error% | Residual |  Time   ")
         println("-------------------")
    end
    for i = 1:options.maxiter
        tic();
        M = method.stepmethod(prob.A,M,options,i);
        timeaccum= timeaccum +  toq(); # Keeps track of time accumulated at every iteration
    
        if(mod(i,options.skip_error_calculation)==0 )
             if(options.exacterror)
                    errors= [ errors vecnorm(prob.Apseudo-M)/initial_error];
             end
             residuals= [residuals vecnorm(prob.A -prob.A*M*prob.A)/initial_residual];
             times = [ times   timeaccum];
             if(options.printiters)
                ## printing iterations info
                @printf "%3.0d  | %3.2f  |  %3.2f  | %3.4f \n" i 100*errors[end] 100*residuals[end] times[end] ;
             end
             if(errors[end] < options.tol)
                fail ="tol-reached"; iterations =i;
             break;
             end
            if(~isempty(options.restol))
                if(residuals[end] < options.restol)
                    fail ="restol-reached"; iterations =i;
                 break;
            end                
            end
            if(isnan(sum(M)) || isnan(errors[end]) || errors[end] >1000  )  
                 fail = "nan";  iterations = i;
             return;
            end
        end 
     if(timeaccum >options.max_time )
         fail ="times_up";  iterations = i;
         if(options.exacterror)
            errors= [ errors vecnorm(prob.Apseudo-M)/initial_error];
         end
         residuals= [residuals vecnorm(prob.A -prob.A*M*prob.A)/initial_residual];
         times = [ times   timeaccum];
         if(options.printiters)
           @printf "%3.0d  | %3.2f  |  %3.2f  | %3.4f \n" i 100*errors[end] 100*residuals[end] times[end] ;
         end
         break;
     end

    end
    if(iterations==0) 
        fail ="max_iter"; 
        iterations=options.maxiter;
    end 
    output = Output(iterations,method.flopsperiter, times, errors, residuals,method.name,fail); 

return output
    
end
