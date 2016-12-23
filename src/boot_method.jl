function boot_method(method_name::AbstractString, prob,options::MyOptions)
@match method_name begin
    "NewtonSchulz"  => method = boot_NewtonSchulz(prob,options);
        "SATAX"     => method = boot_SATAX(prob,options);
            "SAXAS" => method = boot_SAXAS(prob,options);
"NewtonSchulz_warm" => method = boot_NewtonSchulz_warm(prob,options);
                  _ => method = "METHOD DOES NOT EXIST"
    end
    return method;
end

