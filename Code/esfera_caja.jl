#tipos con mayusculas y variables minusculas
#Doctor Watson, jld2

mutable struct esfera{D}#T
    x::SVector{D,Float64}
    v::SVector{D,Float64}
    r::Float64 
    m::Float64 
    c::Float64 
    t::Float64
#     handleT::Int64 podria servir con mutablebinaryminheap
#     handleC::Int64
    tipo::Int64
end

function contarGrCh(Arreglo::Vector{esfera{2}})
    #usar funcion cont
    N = length(Arreglo)
    contador = 0
    contador = count(x->(x.tipo == 1), Arreglo)
    contador, N - contador
end

"""
Documentacion aqui
"""
function leerCaja(nombreDir::String)
    dir = filter(x->occursin(nombreDir,x), readdir(pwd()))[1]
    dfs  = []
    for nombre in ["Datos", "Esferas", "Caja"]
        nombreCompleto = filter(x->occursin(nombre,x), readdir(dir))
        df = CSV.File(dir*"/"*nombreCompleto[1]) |> DataFrame ;
        push!(dfs, df)
    end
    dfs
end
    
function leerEventos(nombreDir::String)
    dir = filter(x->occursin(nombreDir,x), readdir(pwd()))[1]
    nombreCompleto = filter(x->occursin("Eventos",x), readdir(dir))
    df = CSV.File(dir*"/"*nombreCompleto[1]) |> DataFrame ;    
end

function velocidad(D::Int64, vmax::Float64)
    #NO genera aleatorios
    r = rand()*vmax
    ϕ = rand()*2*pi
    # generador de velocidad con distribución circular
    if D == 2
        velocidad = SVector{2,Float64}(r*cos(ϕ), r*sin(ϕ))
    elseif D == 3
        θ = rand()*pi
        velocidad = SVector{3,Float64}(r*sin(θ)*cos(ϕ), r*sin(θ)*sin(ϕ), r*cos(θ))
    end 
end

function velocidadAvg(Esferas::Array{esfera{D}, 1}) where D
    velAvg = zeros(D)
    N = length(Esferas)
    for n1 in 1:N
        velAvg += Esferas[n1].v
    end
    velAvg / N
end

mutable struct vecinosDinamicos
    vecino::Int64
    tiempo::Float64
    ncols::Int64
    tipo::Int64
end

import Base: isless
function isless(a::vecinosDinamicos, b::vecinosDinamicos)
    return a.tiempo < b.tiempo
end

insert!(v::Vector, x) = (splice!(v, searchsorted(v,x), [x]); v)

mutable struct caja{D}
    N::Int64 #numero de esferas
    Ntipos::Int64 #numero de especies
    D::Int64 #dimesion
    esferas::Array{esfera{D}, 1}
    tiempoG::Float64
    propGCrn::Array{Float64,1}
    crecimientos::Array{Float64, 1}
    vecinos::Array{Array{vecinosDinamicos,1},1}
    numeroVecinos::Array{Int64,1}
    numeroColisiones::Array{Int64,1}
    tipoDeVecinos::Array{Int64,2}
    transferenciasPeriodicas::Array{Array{Float64,1},1}
    # potencialSuave::potencial
    #phi::Float64
    caja::Array{Float64,1} #caja cuadrada del espacio de contencion
    L::Float64
end
    
function caja(N::Int, ntipos::Int, D::Int, cr::Float64; propGCrn::Array{Float64, 1} = [1.4, 2.0], vmax::Float64 = 5.0, cajon::Array{Float64,1} = [0.0,1.0])
    numeroVecinos = zeros(Int64, N)
    tipoDeVecinos = zeros(Int64, N, 2)
    numeroColisiones = zeros(Int64, N)
    L = abs(cajon[2] - cajon[1])  # tamaño del dominio
    Esferas = Array{esfera{D}, 1}(undef, N)
    vecinos = Array{Array{vecinosDinamicos,1},1}(undef, N)
    transferenciasPeriodicas = Array{Array{Int64,1},1}(undef, N)
#         transferenciasPeriodicas = zeros(N, D)
    crecimientos = [cr]
    if ntipos == 2  
        P = 1 / (1 + propGCrn[2])
        println("P = $P")
    else
        propGCrn[1] = 1.0
    end
    for n1 in 1:N
        x = SVector{D,Float64}(L*rand(D))
        v = velocidad(D, vmax)
        if ntipos == 1
            Esferas[n1] = esfera{D}(x, v, 0.0, 1.0, cr, 0.0, 0)
        elseif ntipos == 2
            if rand() < P #grandes
                Esferas[n1] = esfera{D}(x, v, 0.0, (propGCrn[1]^2)*1.0, propGCrn[1]*cr, 0.0, 2)
            else #chicas
                Esferas[n1] = esfera{D}(x, v, 0.0, 1.0, cr, 0.0, 1)
            end
        else
            error("Solo se aceptan sistemas con 1 o 2 crecimientos")
        end
        
        transferenciasPeriodicas[n1] = zeros(Int, D)
        vecinos[n1] = []
    end

    # potencialSuave = potencial(0.0, 0.0)

    AvgVel = velocidadAvg(Esferas)
    for n1 in 1:N
        Esferas[n1].v -= AvgVel
    end
    caja{D}(N, ntipos, D, Esferas, 0.0, propGCrn, crecimientos, vecinos, numeroVecinos, numeroColisiones, tipoDeVecinos, transferenciasPeriodicas, cajon, L)#, potencialSuave, cajon, L)
end

# function esferasSuaves(Caja::caja, lambda::Float64, epsilon::Float64)
#     Caja.potencialSuave = potencial(lambda, epsilon)
# end

#     function caja(nombreDir::String)
#         dfs = leerCaja(nombreDir)
#         N = dfs[1].N[1]
#         numeroVecinos = zeros(Int64, N)
#         numeroColisiones = zeros(Int64, N)
#         transferenciasPeriodicas = zeros(Int64, N, 2)
#         Esferas = Array{esfera, 1}()
#         vecinos = Array{Array{Int64,1},1}(undef, N)
#         Ntipos = dfs[1].Ntipos[1]
#         D = dfs[1].D[1]
#         tiempoG = dfs[1].tiempoG[1]
#         propGrCh = dfs[1].propGrCh[1]
#         cajon = [dfs[1].cajon[1], dfs[1].cajon[2]]
#         L = dfs[1].L[1]
#         crecimientos = [dfs[1].crecimiento[1]]
#         for n1 in 1:200
# #             @show dfs[2][n1,:]
#             pos = SVector{2,Float64}(eval(Meta.parse(dfs[2][n1,:].posicion)))
#             vel = SVector{2,Float64}(eval(Meta.parse(dfs[2][n1,:].velocidad)))
#             s = esfera(pos, vel, dfs[2][n1,:].radio, dfs[2][n1,:].masa, dfs[2][n1,:].crecimiento, 
#                 dfs[2][n1,:].tiempo, dfs[2][n1,:].tipo)
#             push!(Esferas, s)
#             transferenciasPeriodicas[n1,:] = [dfs[3][n1,:].transferenciasPeriodicasx, dfs[3][n1,:].transferenciasPeriodicasy]
#             numeroVecinos[n1] = dfs[3][n1,:].numeroVecinos
#             numeroColisiones[n1] = dfs[3][n1,:].numeroColisiones
#             vecinos[n1] = eval(Meta.parse(dfs[3][n1,:].vecinos))
#         end
#         new(N, Ntipos, D, Esferas, tiempoG, propGrCh, crecimientos, vecinos, numeroVecinos, numeroColisiones, transferenciasPeriodicas, cajon, L)
#     end
#     function caja(N::Int64, radio::Float64, posiciones::Array{Float64,2})
#         numeroVecinos = zeros(Int64, N)
#         numeroColisiones = zeros(Int64, N)
#         transferenciasPeriodicas = zeros(Int64, N, 2)
#         cajon = [0.0, 1.0]
#         L = abs(cajon[2] - cajon[1])
#         Esferas = Array{esfera, 1}()
#         vecinos = Array{Array{Int64,1},1}(undef, N)
#         vmax = 5.0
#         for n1 in 1:N
#             posicion = SVector{2,Float64}(posiciones[n1,:])
#             #Velocidades
#             r = rand()*vmax
#             ϕ = rand()*2*pi
#             # generador de velocidad con distribución circular
#             velocidad = SVector{2,Float64}(r*cos(ϕ), r*sin(ϕ))
#             push!(Esferas, esfera(posicion, velocidad, radio, 1.0, 0.0, 0.0, 0))
#             vecinos[n1] = []
#         end
#         new(N, 1, 2, Esferas, 0.0, 1.0, [0.0], vecinos, numeroVecinos, numeroColisiones, transferenciasPeriodicas, cajon, L)
#     end
# end

function correcionPeriodica(X, L::Float64)
    if X > L / 2
        return -L
    elseif X < -L / 2
        return L
    else
        return 0.0
    end
end

function distanciaPorComponente(X, Caja::caja{D}) where D
    X += correcionPeriodica.(X, Caja.L)
end


