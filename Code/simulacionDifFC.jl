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
include(dir*"nuevasMediciones.jl")
# include(dir*"entropy.jl")

ENV["GKSwstype"] = "nul"

normalize(X) = X / norm(X)         

@time medirnSoloDifusion(NUMERO, TIPO, DIMENSION, YAMETE, NANI) # DIF2DT2P
