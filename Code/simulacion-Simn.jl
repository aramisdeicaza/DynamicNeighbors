using LinearAlgebra, StaticArrays
using CSV, DataFrames, Dates
using Plots#, Plots.PlotMeasures, ColorSchemes, LaTeXStrings#, Makie
using LsqFit, SphericalHarmonics
using DataStructures
default(size = (900, 600), titlefont = (10, "times"), guidefont = (13, "times"), tickfont = (10, "times"))
#     legendfont = (12, "times"))#, yguidefontrotation = true) 
dir = "/home/aramis/ICN-PythiaTest/Codigo/"

include(dir*"esfera_caja.jl")
# include(dir*"triangula_fcc.jl")
include(dir*"celdas.jl")
include(dir*"eventos_tiempo.jl")
include(dir*"visual.jl")
include(dir*"generador_iterador.jl")
include(dir*"mediciones.jl")
# include(dir*"entropy.jl")

ENV["GKSwstype"] = "nul"

normalize(X) = X / norm(X)         

prop=[RADIO, PROPORCION]
dt = 1e-10 

@time medirSimetriasAtascadas(10044, 2, 2, 8e-2, dt, propGCrn=prop) # SIMULACION