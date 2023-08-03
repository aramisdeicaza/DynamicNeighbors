using LinearAlgebra, StaticArrays
using CSV, DataFrames, Dates
using Plots#, Plots.PlotMeasures, ColorSchemes, LaTeXStrings#, Makie
using LsqFit, SphericalHarmonics
using DataStructures
default(size = (900, 600), titlefont = (10, "times"), guidefont = (13, "times"), tickfont = (10, "times"))
#     legendfont = (12, "times"))#, yguidefontrotation = true) 

dir = "/storage/alice/Edgar_alice/ICN-PythiaTest/Codigo/"

include(dir*"esfera_caja.jl")
# include(dir*"triangula_fcc.jl")
include(dir*"celdas.jl")
include(dir*"eventos_tiempo.jl")
include(dir*"visual.jl")
include(dir*"generador_iterador.jl")
include(dir*"mediciones.jl")
# include(dir*"entropy.jl")

ENV["GKSwstype"] = "nul"
# N hexatic 256^2 = 65_536 o 512^2 = 262_144, 1024^2
# Crecimientos2d 2k-10k 8e-3 - Crecimientos3d 2k rapido 1e-1 lento 2e-3

# if NUMERO == 2048
#     tdifusion = 19.0
# end

#if DIMENSION == 2
#    gr()
#elseif DIMENSION == 3
#    pyplot()
#end

normalize(X) = X / norm(X)         

# @time medirSoloDifusion(NUMERO, TIPO, DIMENSION, YAMETE, NANI, tdifusion) # DIFUSION
@time medirTodo(NUMERO, TIPO, DIMENSION, YAMETE, NANI) # SIMULACION


# Falta script para esto
# @time medirSoloDifusion(NUMERO, 1, 2, YAMETE, NANI, propGrCh = NANDATO) # DIF2DT2P
# @time medirTodo(NUMERO, 2, 2, YAMETE, NANI, propGrCh = NANDATO) # SIM2DT2P
