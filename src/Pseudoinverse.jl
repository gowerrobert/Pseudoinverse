module Pseudoinverse
using StatsBase
using Match
using Plots
#using Winston

# The types exported
export MyOptions, Prob, Method, Output
# The exported functions
export pseudoinvert_matrix, uniform_mat_rank, uniform_sym_rank
export plot_outputs_Plots

type MyOptions
    n::Int
    m::Int
    sketchsize::Int
    sketch::AbstractString
    tol::Float64
    restol::Float64
    maxiter::Int
    skip_error_calculation::Int
    max_time::Float64
    printiters::Bool 
    exacterror::Bool
    M0type::AbstractString
end


type Prob
    A::Array{Float64}
    b::Array{Float64}
    xsol::Array{Float64}
    Apseudo::Array{Float64}
    name::AbstractString
end

type Method
    flopsperiter::Int
    name::AbstractString
    M0::Array{Float64}
    stepmethod::Function
end

type Output
    iterations::Int
    flopsperiter::Int
    times::Array{Float64}
    errors::Array{Float64}
    residuals::Array{Float64}
    name::AbstractString
    fail::AbstractString
end
#Including method wrappers 
include("pseudoinvert_matrix.jl")
include("boot_method.jl")
#Including test and problem generating functions
include("uniform_mat_rank.jl") 
include("uniform_sym_rank.jl") 
#Including iterative methods for calculating pseudoinverse
include("NewtonSchulz.jl") # NewtonSchulz
include("NewtonSchulz_warm.jl")  #NewtonSchulz with randomized warm starting 
include("SATAX.jl")
include("SAXAS.jl")
#Including utilities, plotting, data analysis
include("plot_outputs_Plots.jl")


end