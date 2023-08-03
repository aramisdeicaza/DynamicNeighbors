using Base: Float64
# include("esfera_caja.jl")
# include("triangula_fcc.jl")
# include("celdas.jl")
# include("eventos_tiempo.jl")
# include("visual.jl")
# include("generador_iterador.jl")

function packingFraction(Caja::caja{D}) where D
    phi = 0.0
    for n1 in 1:Caja.N
        phi += (Caja.esferas[n1].r + (Caja.tiempoG - Caja.esferas[n1].t)*Caja.esferas[n1].c )^D
    end
    phi *= π
    if D == 3
        phi *= (4/3)
    end
    phi
end

function g(Caja::caja{D}, rs::Array{Float64,1}, radio::Float64; maximo::Int64 = 20) where D
    g = zeros(length(rs))
    dr = rs[2] - rs[1]
    for n1 in 1:Caja.N-1
        for n2 in (n1+1):Caja.N
            X = Caja.esferas[n1].x - Caja.esferas[n2].x
            X = distanciaPorComponente(X, Caja)
            if norm(X) < maximo*radio
                j::Int64 = ceil(norm(X)/dr)
                g[j] += 2
            end
        end
    end
    for k in 1:length(rs)
        if D == 2 
            g[k] = g[k] / (Caja.N * 2π * rs[k] )
        elseif D == 3
            g[k] = g[k] / (Caja.N * 4π * (rs[k])^2 )
        end
    end
    g[1] = 0.0 
    g
end

function gbinaria(Caja::caja{2}, rs::Array{Float64,1}, radios::Array{Float64, 1}; maximo::Int64 = 20)
    g = zeros(length(rs))
    g11 = zeros(length(rs))
    g12 = zeros(length(rs))
    g22 = zeros(length(rs))
    dr = rs[2] - rs[1]
    for n1 in 1:Caja.N-1
        for n2 in (n1+1):Caja.N
            X = Caja.esferas[n1].x - Caja.esferas[n2].x
            X = distanciaPorComponente(X, Caja)
            if norm(X) < maximo*radios[2]
                j::Int64 = ceil(norm(X)/dr)
                g[j] += 2
                if Caja.esferas[n1].tipo == 1 && Caja.esferas[n2].tipo == 1
                    g11[j] +=2
                elseif Caja.esferas[n1].tipo == 1 && Caja.esferas[n2].tipo == 2
                    g12[j] += 2
                elseif Caja.esferas[n1].tipo == 2 && Caja.esferas[n2].tipo == 1
                    g12[j] += 2
                elseif Caja.esferas[n1].tipo == 2 && Caja.esferas[n2].tipo == 2
                    g22[j] += 2
                end
            end
        end
    end
    for k in 1:length(rs)
        g[k] = g[k] / (Caja.N * 2π * rs[k] )
        g11[k] = g11[k] / (Caja.N * 2π * rs[k] )
        g12[k] = g12[k] / (Caja.N * 2π * rs[k] )
        g22[k] = g22[k] / (Caja.N * 2π * rs[k] )
    end
    g[1], g11[1], g12[1], g22[1] = (0.0, 0.0, 0.0, 0.0)
    g, g11, g12, g22
end       
                                                              
function distanciasPeriodicas(n1::Int64, n2::Int64, Caja::caja{D}) where D
    x1 = Caja.esferas[n1].x
    x2 = Caja.esferas[n2].x
    Xs = []
    CIP = [] # correccionImagenesPeriodicas
    if D == 2
        for h in -1:1, l in -1:1
            c = [h,l]
            push!(CIP, c)
        end
    elseif D == 3
         for h in -1:1, l in -1:1, g in -1:1
            c = [h,l,g]
            push!(CIP, c)
        end   
    end
    for xp in CIP
        push!(Xs, norm(x2+xp-x1))
    end
    Xs
end

function gdp(Caja::caja{D}) where D
    dr = 1 / (2*Caja.N)
    rs = collect(0:dr:5.0)
    g = zeros(length(rs))
    
    for n1 in 1:Caja.N-1
        for n2 in (n1+1):Caja.N
            Xs = distanciasPeriodicas(n1, n2, Caja)
            for X in Xs
                j::Int64 = ceil(norm(X)/dr)
                g[j] += 2
            end
        end
    end
    for k in 1:length(rs)
        if D == 2 
            g[k] = g[k] / (9*Caja.N * 2π * rs[k] )
        elseif D == 3
            g[k] = g[k] / (9*Caja.N * 4π * (rs[k])^2 )
        end
    end
    g[1] = 0.0
    g, rs ./ Caja.esferas[1].r
end

function gbinariaP(Caja::caja{2})
    dr = 1 / (2*Caja.N)
    rs = collect(0:dr:5.0)
    
    g = zeros(length(rs))
    g11 = zeros(length(rs))
    g12 = zeros(length(rs))
    g22 = zeros(length(rs))
    
    for n1 in 1:Caja.N-1
        for n2 in (n1+1):Caja.N
            
            Xs = distanciasPeriodicas(n1, n2, Caja)
            for X in Xs
                j::Int64 = ceil(norm(X)/dr)
                g[j] += 2
                if Caja.esferas[n1].tipo == 1 && Caja.esferas[n2].tipo == 1
                    g11[j] +=2
                elseif Caja.esferas[n1].tipo == 1 && Caja.esferas[n2].tipo == 2
                    g12[j] += 2
                elseif Caja.esferas[n1].tipo == 2 && Caja.esferas[n2].tipo == 1
                    g12[j] += 2
                elseif Caja.esferas[n1].tipo == 2 && Caja.esferas[n2].tipo == 2
                    g22[j] += 2
                end
            end
        end
    end
    for k in 1:length(rs)
        g[k] = g[k] / (9*Caja.N * 2π * rs[k] )
        g11[k] = g11[k] / (9*Caja.N * 2π * rs[k] )
        g12[k] = g12[k] / (9*Caja.N * 2π * rs[k] )
        g22[k] = g22[k] / (9*Caja.N * 2π * rs[k] )
    end
    g[1], g11[1], g12[1], g22[1] = (0.0, 0.0, 0.0, 0.0)
                
    radioMin = rminmax(Caja)
    g, g11, g12, g22, rs ./ radioMin[1]
end

function medirDistribucionRadialP(Caja::caja{D}, Celdas::celdas{D}, Eventos::eventos, phi::Float64, crecimiento::Float64; modNombre::String="", max::Int64 = 20) where D
    DatosNombreArchivo = datosNombreArchivo(Caja.N, Caja.D, Caja.Ntipos, phi, crecimiento)
    NombreArchivoGrS = "GrPS"*modNombre*DatosNombreArchivo*".csv" 
    NombreArchivoGr = "GrP"*modNombre*DatosNombreArchivo*".csv" 
                                                                            
    Ca, Ce, Ev = (deepcopy(Caja), deepcopy(Celdas), deepcopy(Eventos))

    Gt, rs = gdp(Ca) 
    CSV.write(NombreArchivoGrS, DataFrame(distancia = rs, g = Gt))

    prints = 1
    n::Int64 = cld(Ca.N, 10)
    for i in 1:999
        if prints / 5 <= i / 999
            println("Avance de gP(r) a $prints / 5 ")
            prints += 1
        end
        Caja, Celdas, Eventos = evolucionarNumeroPasos(Ca, Ce, Ev, n)
        Caja = mover_esferas(Ca)
        G, rss = gdp(Ca)
        Gt += G
    end
    Gt /= 1000
    
    CSV.write(NombreArchivoGr, DataFrame(distancia = rs, g = Gt))
    
    plot(rs, Gt)
    plot!(xlim=(0,50), label = "phi = $phi", xlabel = "r", ylabel = "gp(r)")
    savefig("gr"*modNombre*"-"*DatosNombreArchivo*".png")
end

function medirDistribucionRadialBinarioP(Caja::caja{2}, Celdas::celdas{2}, Eventos::eventos, phi::Float64, crecimiento::Float64;  modNombre::String="", max::Int64 = 20)
    DatosNombreArchivo = datosNombreArchivo(Caja.N, Caja.D, Caja.Ntipos, phi, crecimiento)
    NombreArchivoGrS = "GrPS"*modNombre*DatosNombreArchivo*".csv" 
    NombreArchivoGr = "GrP"*modNombre*DatosNombreArchivo*".csv" 
                                                                                        
   
    Ca, Ce, Ev = (deepcopy(Caja), deepcopy(Celdas), deepcopy(Eventos))
                                                                            
    G, G11, G12, G22, rs = gbinariaP(Ca)
    CSV.write(NombreArchivoGrS, DataFrame(distancia = rs, g = G, g11 = G11, g12 = G12, g22 = G22))
       
    prints = 1
    n::Int64 = cld(Ca.N,  10)
    for i in 1:999
        if prints / 5 <= i / 999
            println("Avance de gP(r) a $prints / 5 ")
            prints += 1
        end
        Ca, Ce, Ev = evolucionarNumeroPasos(Ca, Ce, Ev, n)
        Ca = mover_esferas(Ca)
        
        Gs, G11s, G12s, G22s, rss = gbinariaP(Ca)
        G += Gs
        G11 += G11s
        G12 += G12s
        G22 += G22s
    end
    G /= 1000
    G11 /= 1000
    G12 /= 1000
    G22 /= 1000
    
    CSV.write(NombreArchivoGr, DataFrame(distancia = rs, g = G, g11 = G11, g12 = G12, g22 = G22))
    
    plot(rs, G, label = "g")
    plot!(rs, G11, label = "g11")
    plot!(rs, G12, label = "g12")
    plot!(rs, G22, label = "g22")
    plot!(xlim=(0,50), title = "phi = $phi", xlabel = "r", ylabel = "gp(r)")
    savefig("gr"*modNombre*"-"*DatosNombreArchivo*".png")
end
                                    
# struct vecinoDistancia
#     vecino::Int64
#     ri::Float64
# end

# function isless(a::vecinoDistancia, b::vecinoDistancia)
#     return a.ri < b.ri
# end
                                    
function VecinosParametrodeOrden(Caja::caja{D}, r::Float64=3.0) where D #Solo para comprabar que el otro algoritmo funcione
    N = Caja.N
#     VecinosRi = [ BinaryMinHeap{vecinoDistancia}() for i in 1:N]
    Vecinos = [ Int64[] for i in 1:N]
    Distancias = [ Float64[] for i in 1:N]
    
    R = 1 + r
    
    for i in 1:N-1
        for j in i+1:N
            Δp = Caja.esferas[i].x - Caja.esferas[j].x
            Δp = distanciaPorComponente(Δp, Caja) 
            ri = norm(Δp)
                
            if ri < (Caja.esferas[i].r*R + Caja.esferas[j].r)
                                                    
                push!(Vecinos[i], j)
                push!(Vecinos[j], j)
                push!(Distancias[i], ri)
                push!(Distancias[j], ri)
            end
        end
        perm = sortperm(Distancias[i])
        Vecinos[i] = Vecinos[i][perm]
#         sort!(Vecinos[i], by=x->norm(distanciaPorComponente((Cajas.esferas[x].x-Caja.esferas[i].x), Caja)))
    end
#     sort!(Vecinos[N], by=x->norm(distanciaPorComponente((Cajas.esferas[x].x-Caja.esferas[N].x), Caja)))
    perm = sortperm(Distancias[N])
    Vecinos[N] = Vecinos[N][perm]
                                        
    nv = 0
    if D == 2
        nv = 6
    elseif D == 3
        nv = 12
    end
                                            
    VPO = [ Int64[] for i in 1:N]
    for i in 1:N
        if length(Vecinos[i]) >= nv
            VPO[i] = Vecinos[i][1:nv]
        else
            VPO[i] = Vecinos[i]
        end
    end
    VPO
end

function angulo(n1::Int64, n2::Int64, Caja::caja{2})
    X = Caja.esferas[n1].x - Caja.esferas[n2].x
    X = distanciaPorComponente(X, Caja)
    θ = atan(X[2], X[1])
end

function psi6(Caja::caja{2}; ncols::Int64=1) 
    PO6 = 0
    for n1 in 1:Caja.N
        psi = 0
        for n2vd in Caja.vecinos[n1]
            if n2vd.ncols >= ncols
                psi += exp(im*6*angulo(n1, n2vd.vecino, Caja))
            end
        end
        nv = count(x->(x.ncols >= ncols), Caja.vecinos[n1])
        if nv > 2
            psi /= nv
        else
            psi = 0
        end
        PO6 += norm(psi)
    end
    PO6 /= Caja.N
    PO6 
end

function psi(Caja::caja{2}, orden; ncols::Int64=1) 
    PO = 0
    for n1 in 1:Caja.N
        psi = 0
        for n2vd in Caja.vecinos[n1]
            if n2vd.ncols >= ncols
                psi += exp(im*orden*angulo(n1, n2vd.vecino, Caja))
            end
        end
        nv = count(x->(x.ncols >= ncols), Caja.vecinos[n1])
        if nv > 2
            psi /= nv
        else
            psi = 0
        end
        PO += norm(psi)
    end
    PO /= Caja.N
    PO
end

function meanSquareDisplacement(spheres0::Array{esfera{D}, 1},  Caja::caja{D}) where D
    MSD = 0
    for n1 in 1:Caja.N
        MSD += (norm(Caja.esferas[n1].x + Caja.transferenciasPeriodicas[n1] - spheres0[n1].x))^2
    end
    MSD /= Caja.N
    MSD
end                         
                        
function meanSquareDisplacement(spheres0::Array{esfera{D}, 1}, velAvg, dt, Caja::caja{D}) where D
    MSD = 0
    for n1 in 1:Caja.N
        MSD += (norm(Caja.esferas[n1].x + Caja.transferenciasPeriodicas[n1] - velAvg*dt - spheres0[n1].x))^2
        # MSD += (Caja.esferas[n1].x + Caja.transferenciasPeriodicas[n1] - spheres0[n1].x)⋅(Caja.esferas[n1].x + Caja.transferenciasPeriodicas[n1] - spheres0[n1].x)
    end
    MSD /= Caja.N
    MSD
end 

function polarAngle(n1::Int64, n2::Int64, Caja::caja{3}) #θ ∈ [0, π]
    X = Caja.esferas[n1].x - Caja.esferas[n2].x
    X = distanciaPorComponente(X, Caja)
    θ = atan(sqrt(X[1]^2 + X[2]^2), X[3])
end    

function azimuthAngle(n1::Int64, n2::Int64, Caja::caja{3}) #θ ∈ [-π, π]
    X = Caja.esferas[n1].x - Caja.esferas[n2].x
    X = distanciaPorComponente(X, Caja)
    ϕ = atan(X[2], X[1])
end

function Q6(Caja::caja{3}; ncols::Int64=1) 
    PO = 0
    constant = (4*pi/13)
    for m in -6:6
        q6 = 0
        k = 0
        for n1 in 1:Caja.N
            for n2vd in Caja.vecinos[n1]
                if n2vd.ncols >= ncols
                    k += 1
                    ϕ = azimuthAngle(n1, n2vd.vecino, Caja)
                    θ = polarAngle(n1, n2vd.vecino, Caja)
                    Y = computeYlm(θ, ϕ, lmax = 6) 
                    q6 += Y[(6,m)]
                end
            end
        end
        PO += (norm(q6/k))^2
    end
    PO = sqrt(constant*PO)
    PO / 0.575 #normalize to fcc
end
                                                
function psi6(Caja::caja{2}, vecinosPO::Array{Array{Int64,1},1}) 
    PO6 = 0
    for n1 in 1:Caja.N
        psi = 0
        for n2 in vecinosPO[n1]
            psi += exp(im*6*angulo(n1, n2, Caja))
        end
        nv = length(vecinosPO[n1])
        psi /= 6
        PO6 += norm(psi)
    end
    PO6 /= Caja.N
    PO6 
end
    
function Q6(Caja::caja{3}, vecinosPO::Array{Array{Int64,1},1}) 
    PO = 0
    constant = (4*pi/13)
    for m in -6:6
        q6 = 0
        k = 0
        for n1 in 1:Caja.N
            for n2 in vecinosPO[n1]
                ϕ = azimuthAngle(n1, n2, Caja)
                θ = polarAngle(n1, n2, Caja)
                Y = computeYlm(θ, ϕ, lmax = 6) 
                q6 += Y[(6,m)]
            end
        end
        PO += (norm(q6/12))^2
    end
    PO = sqrt(constant*PO)
    PO / 0.575 #normalize to fcc
end

function rapidez(Caja::caja{D}) where D
    vels = zeros(Caja.N)
    for n1 in 1:Caja.N
        vels[n1] = sqrt(Caja.esferas[n1].v⋅Caja.esferas[n1].v)
    end
    vels
end

function energia(Caja::caja{D}) where D
    Energia  = 0
    for n1 in 1:Caja.N
        Energia += Caja.esferas[n1].m*(Caja.esferas[n1].v⋅Caja.esferas[n1].v) / 2
    end
    Energia
end

function rminmax(Caja::caja{D}) where D
    if Caja.esferas[1].c != 0.0 
        error("Aun en crecimiento")
    end
    rmima = [0.0, 0.0]
    rmin = true
    rmax = true
    if Caja.Ntipos == 1
        return [Caja.esferas[1].r, Caja.esferas[1].r]
    elseif Caja.Ntipos == 2
        for n1 in 1:Caja.N
            if rmin && Caja.esferas[n1].tipo == 1
                rmima[1] = Caja.esferas[n1].r
                rmin = false
            elseif rmax && Caja.esferas[n1].tipo == 2
                rmima[2] = Caja.esferas[n1].r
                rmax = false
            elseif !rmin && !rmax 
                break
            end
        end
    end
    rmima
end

function histograma(Arreglo::Array{Int64,1}, nbins = 100)
    delta = 1#length(Arreglo) / nbins
    
    centros = zeros(nbins)
    conteos = zeros(nbins)
    
    start = 0
    for k in 1:nbins
        stop = start + delta
        conteos[k] = count(i -> start <= i < stop, Arreglo)
        if delta == 1
            centros[k] = start #+ delta/2.
        else
            centros[k] = start + delta/2
        end
        start = stop
   end
        
   conteos, centros 
end

function histogramaMono(delta::Int64, Caja::caja{D}) where D
    nbins::Int = Caja.N/delta
    histograma(Caja.numeroVecinos, nbins)
end


# function histogramaBi(delta::Int64, Caja::caja{D}) where D
#     nbins::Int = Esferas.N/delta
#     conteos = zeros(nbins)
#     conteosChicas = zeros(nbins)
#     conteosGrandes = zeros(nbins)
#     centros = zeros(nbins)
    
#     start = 0
#     for k in 1:nbins
#         stop = start + delta
#         conteos[k] = count(i -> start <= i < stop, Esferas.numeroVecinos)
#         cCh = 0
#         cGr = 0
        
#         rCh = minimum(Esferas.radios)
#         rGr = maximum(Esferas.radios)
#         for i in 1:Esferas.N
#             if start <= Esferas.numeroVecinos[i] < stop && Esferas.radios[i] == rCh
#                 cCh += 1
#             elseif start <= Esferas.numeroVecinos[i] < stop && Esferas.radios[i] == rGr
#                 cGr += 1
#             end
#         end
        
#         conteosChicas[k] = cCh
#         conteosGrandes[k] = cGr
#         conteos[k] == cCh + cGr ? nothing : error("Suma de chicas y grandes no es la total")
        
#         if delta == 1
#             centros[k] = start #+ delta/2.
#         else
#             centros[k] = start + delta/2
#         end
#         start = stop
#    end
        
#    conteos, conteosChicas, conteosGrandes, centros
# end

function gaussDis(Nesferas::Array{Float64,1}, Nvecinos::Array{Float64,1}, promedio::Float64)
    @show length(Nesferas), length(Nvecinos), promedio
    @. model(x, p) = 1 / ( (p[1]*sqrt(2*pi))*exp(-0.5*(x - p[2])^2 / p[1]^2))
    p0 = [1.0, promedio]       # adivinanza inicial
    fit = curve_fit(model, Nvecinos, Nesferas, p0)
    fit.param
end

function medirEnergia(Caja::caja{D}, Celdas::celdas{D}, Eventos::eventos, phi::Float64, crecimiento::Float64; modNombre::String="") where D
    # n es numero de muestras    
    DatosNombreArchivo = datosNombreArchivo(Caja.N, Caja.D, Caja.Ntipos, phi, crecimiento)
    NombreArchivoEnergia = "Energy"*modNombre*DatosNombreArchivo*".csv"
    
    Energia = energia(Caja)
    VelAvg = sum(rapidez(Caja))/ Caja.N
    
    CSV.write(NombreArchivoEnergia, DataFrame(energia = Energia, sqrtE = sqrt(Energia), velAvg = VelAvg))
                
    sqrt(Energia)
end

function medirVariacionEnergia(Caja::caja{D}, Celdas::celdas{D}, Eventos::eventos, phi::Float64, crecimiento::Float64) where D
    # n es numero de muestras
    n = 10_0000
    DatosNombreArchivo = datosNombreArchivo(Caja.N, Caja.D, Caja.Ntipos, phi, crecimiento)
    NombreArchivoEnergia = "Energy"*DatosNombreArchivo*".csv"
    
    Energias = zeros(n)
    Energias[1] = energia(Esferas)
    Tiempos = zeros(n)
    Tiempos[1] = Caja.tiempoG
    prints = 1
    for i in 1:(n)-1
        if prints / 5 < i / (n)
            println("Avance de calculo de Energia a $prints / 5 ")
            prints += 1
        end
        Caja, Celdas, Eventos = evolucionarNumeroPasos(Caja, Celdas, Eventos, CajaN)
        Energias[i+1] = energia(Caja)
        Tiempos[i+1] = Esferas.tiempoG
    end
    
    CSV.write(NombreArchivoEnergia, DataFrames(energia = Energias, tiempo = Tiempos))
    
    plot(Tiempos, Energias)
    plot!(label = "phi = $phi", xlabel = "tiempo", ylabel = "Energia")
    savefig("E-"*DatosNombreArchivo*".png")
end

function medirDistribucionVecinosDinamicos(Caja::caja{D}, Celdas::celdas{D}, Eventos::eventos, phi::Float64, crecimiento::Float64; modNombre::String="") where D
    DatosNombreArchivo = datosNombreArchivo(Caja.N, Caja.D, Caja.Ntipos, phi, crecimiento)
    NombreArchivoHist = "Hist"*modNombre*DatosNombreArchivo*".csv"
    NombreArchivoAvg = "Avg"*modNombre*DatosNombreArchivo*".csv"

    nesferas, nvecinos = histogramaMono(1, Caja)
    AvgVecinos = sum(Caja.numeroVecinos)/Caja.N
    AvgColisiones = sum(Caja.numeroColisiones)/Caja.N
    FitGauss = 0.0 #gaussDis(nesferas, nvecinos, AvgVecinos)
    #Genera archivo de histograma
    CSV.write(NombreArchivoHist, DataFrame(nvecinos = nvecinos, nesferas = nesferas))
    
    #Genera archivos con promedio de vecinos y ajuste de gauss
    if Caja.Ntipos == 1
        AvgVecinos1, AvgVecinos2 = (0.0, 0.0)
    elseif Caja.Ntipos == 2
        AvgVecinos1 = sum(count.(x->(x.tipo == 1), Caja.vecinos)) / Caja.N
        AvgVecinos2 = AvgVecinos - AvgVecinos1
    end
    CSV.write(NombreArchivoAvg, DataFrame(AvgVecinos = AvgVecinos, AvgColisiones = AvgColisiones, FitGauss = FitGauss, AvgVecinos1 = AvgVecinos1, AvgVecinos2 = AvgVecinos2))
end

function medirPsi6t(Caja::caja{D}, Celdas::celdas{D}, Eventos::eventos, phi::Float64, crecimiento::Float64; modNombre::String="") where D
    DatosNombreArchivo = datosNombreArchivo(Caja.N, Caja.D, Caja.Ntipos, phi, crecimiento)
    NombreArchivoPO6T = "PO6T"*modNombre*DatosNombreArchivo*".csv"
    # Caja = reiniciar(Caja)

    n = 1000
    pos = zeros(n)
    tiempos = LinRange(0, 1, n)
    for i in 1:n
        Caja, Celdas, Eventos = evolucionardt(Caja, Celdas, Eventos, 0.001)
        Caja = mover_esferas(Caja)
        PO = psi6(Caja)
        pos[i] = PO
    end
    
    CSV.write(NombreArchivoPO6T, DataFrame(psi6 = pos, tiempo = tiempos ))
    pos, tiempo
end

function medirParametroOrdenOrientacional(Caja::caja{D}, Celdas::celdas{D}, Eventos::eventos, phi::Float64, crecimiento::Float64; modNombre::String="") where D
    DatosNombreArchivo = datosNombreArchivo(Caja.N, Caja.D, Caja.Ntipos, phi, crecimiento)
    NombreArchivoPO = "PO"*modNombre*DatosNombreArchivo*".csv"

    if D == 2
        PO = psi(Caja, 6), psi(Caja, 12)
        CSV.write(NombreArchivoPO, DataFrame(psi6 = PO[1], psi12=PO[2]))
        
        NombreArchivoPOX = "POX"*modNombre*DatosNombreArchivo*".csv"
        O = 0:0.01:15
        POX = zeros(length(O))
        for i in 1:length(O)
            POX[i] = psi(Caja, O[i])
        end
        CSV.write(NombreArchivoPOX, DataFrame(orden = O, psix = POX))
        plot(title = "Con $Caja.ntipos radios, $Caja.N discos",
            xlabel = "Orden",  ylabel = "psi x", xtick=0:1:15, label = "phi = $phi")
        plot!(O, POX, label = "Todos")
        savefig("pox-"*modNombre*DatosNombreArchivo*".png")
    elseif D == 3 
        PO = Q6(Caja)
        CSV.write(NombreArchivoPO, DataFrame(phi = phi, Q6 = PO))
    end
end

function medirDistribucionRadial(Caja::caja{D}, Celdas::celdas{D}, Eventos::eventos, phi::Float64, crecimiento::Float64; modNombre::String="", max::Int64 = 20) where D
    DatosNombreArchivo = datosNombreArchivo(Caja.N, Caja.D, Caja.Ntipos, phi, crecimiento)
    NombreArchivoGrS = "GrS"*modNombre*DatosNombreArchivo*".csv" 
    NombreArchivoGr = "Gr"*modNombre*DatosNombreArchivo*".csv" 

    radio = Caja.esferas[1].r

    dr = 1 / (2*Caja.N)
    rs = collect(0:dr:max*radio)
    G = g(Caja, rs, radio, maximo = max) 
    CSV.write(NombreArchivoGrS, DataFrame(distancia = rs./radio, gr = G))

    prints = 1
    n::Int64 = cld(Caja.N, 10)
    for i in 1:999
        if prints / 5 <= i / 999
            println("Avance de g(r) a $prints / 5 ")
            prints += 1
        end
        Caja, Celdas, Eventos = evolucionarNumeroPasos(Caja, Celdas, Eventos, n)
        Caja = mover_esferas(Caja)
        G += g(Caja, rs, radio, maximo = max)
    end
    G /= 1000
    rs = rs ./ radio
    
    CSV.write(NombreArchivoGr, DataFrame(distancia = rs, gr = G))
    
    plot(rs, G)
    plot!(xlim=(0,10), label = "phi = $phi", xlabel = "r", ylabel = "g(r)")
    savefig("gr"*modNombre*"-"*DatosNombreArchivo*".png")
end

function medirDistribucionRadialBinario(Caja::caja{2}, Celdas::celdas{2}, Eventos::eventos, phi::Float64, crecimiento::Float64;  modNombre::String="", max::Int64 = 20)
    DatosNombreArchivo = datosNombreArchivo(Caja.N, Caja.D, Caja.Ntipos, phi, crecimiento)
    NombreArchivoGrS = "GrS"*modNombre*DatosNombreArchivo*".csv" 
    NombreArchivoGr = "Gr"*modNombre*DatosNombreArchivo*".csv" 
    
    dr = 1 / (2*Caja.N)
    radios = rminmax(Caja)

    rs = collect(0:dr:max*radios[2])
    G, G11, G12, G22 = gbinaria(Caja, rs, radios, maximo=max) 
    rs1 = rs ./ radios[1]
    rs2 = rs ./ radios[2]
    CSV.write(NombreArchivoGrS, DataFrame(distancia = rs1, distancia2 = rs2, g = G, g11 = G11, g12 = G12, g22 = G22))
       
    prints = 1
    n::Int64 = cld(Caja.N,  10)
    for i in 1:999
        if prints / 5 <= i / 999
            println("Avance de g(r) a $prints / 5 ")
            prints += 1
        end
        Caja, Celdas, Eventos = evolucionarNumeroPasos(Caja, Celdas, Eventos, n)
        Caja = mover_esferas(Caja)
        
        Gs, G11s, G12s, G22s = gbinaria(Caja, rs, radios, maximo=max)
        G += Gs
        G11 += G11s
        G12 += G12s
        G22 += G22s
    end
    G /= 1000
    G11 /= 1000
    G12 /= 1000
    G22 /= 1000
    rs1 = rs ./ radios[1]
    rs2 = rs ./ radios[2]
    
    CSV.write(NombreArchivoGr, DataFrame(distancia = rs1, distancia2 = rs2, g = G, g11 = G11, g12 = G12, g22 = G22))
    
    plot(rs1, G, label = "g")
    plot!(rs1, G11, label = "g11")
    plot!(rs1, G12, label = "g12")
    plot!(rs1, G22, label = "g22")
    plot!(xlim=(0,10), title = "phi = $phi", xlabel = "r", ylabel = "g(r)")
    savefig("gr"*modNombre*"-"*DatosNombreArchivo*".png")
end

function medirDifusion(Caja::caja{D}, Celdas::celdas{D}, Eventos::eventos, crecimiento::Float64, phi::Float64, tdif::Real) where D
    DatosNombreArchivo = datosNombreArchivo(Caja.N, Caja.D, Caja.Ntipos, phi, crecimiento)
    NombreArchivoD = "Dif"*DatosNombreArchivo*".csv"
    NombreArchivoPs = "P12-"*DatosNombreArchivo*".csv"

    esferasIniciales = deepcopy(Caja.esferas)
    tiempoInicial = Caja.tiempoG
    tiempoFinal = tiempoInicial + tdif
    
    MSDS = Float64[]
    MSDVS = Float64[]
    tiempos = Float64[]
    colisiones = Float64[]
    
    P1 = []
    T1 = []
    P2 = []
    T2 = []
#     velAvg = velocidadAvg(Caja.esferas)
    prints = 1
    while Caja.tiempoG < tiempoFinal
        if prints*tdif/5 < (Caja.tiempoG - tiempoInicial)
            println("Avance de difusion a $prints / 5")
            prints += 1
        end
        Caja, Celdas, Eventos = evolucionardt(Caja, Celdas, Eventos, 0.001, addvecinos=false)
        Caja = mover_esferas(Caja)
        
        dt = Caja.tiempoG - tiempoInicial
        MSD = meanSquareDisplacement(esferasIniciales, Caja)
            
        push!(MSDS, MSD)
        push!(MSDVS, MSD / (sum(rapidez(Caja))^2 / Caja.N) )
        push!(tiempos, dt)
        push!(colisiones, sum(Caja.numeroColisiones)/Caja.N)
                                                                                    
        push!(P1, Caja.esferas[1].x)                                                
        push!(P2, Caja.esferas[2].x)
        push!(T1, Caja.transferenciasPeriodicas[1])                                                
        push!(T2, Caja.transferenciasPeriodicas[2])
    end
    CSV.write(NombreArchivoD, DataFrame(MSD = MSDS, MSDV = MSDVS, tiempo = tiempos, colisiones = colisiones ))
    CSV.write(NombreArchivoPs, DataFrame(P1 = P1, T1 = T1, P2 = P2, T2 = T2, tiempo = tiempos))
end

function medirLineasDeTiempo(Caja::caja{D}, Celdas::celdas{D}, Eventos::eventos, phi::Float64, crecimiento::Float64, tldt::Float64) where D
    DatosNombreArchivo = datosNombreArchivo(Caja.N, Caja.D, Caja.Ntipos, phi, crecimiento)
    NombreArchivoLT = "LT"*DatosNombreArchivo*".csv"

    AvgVecinos = Float64[]
    AvgColisiones = Float64[]
    tiempos = Float64[]

    tiempoInicial = Caja.tiempoG
    tiempoFinal = tiempoInicial + tldt
    
    n::Int64 = cld(Caja.N,  10)
    prints = 1
    while Caja.tiempoG < tiempoFinal
        if prints*tldt/5 < (Caja.tiempoG - tiempoInicial)
            println("Avance de difusion a $prints / 5")
            prints += 1
        end
        Caja, Celdas, Eventos = evolucionarNumeroPasos(Caja, Celdas, Eventos, n)
        push!(AvgVecinos, sum(Caja.numeroVecinos)/Caja.N)
        push!(AvgColisiones, sum(Caja.numeroColisiones)/Caja.N)
        push!(tiempos, Caja.tiempoG)
    end
    CSV.write(NombreArchivoLT, DataFrame(AvgVecino = AvgVecinos, tiempo = tiempos))
end

function medirTodoTriFCC()
    DatosNombreArchivo = datosNombreArchivo(N, dimension, tipo, phi, crecimiento)
    

end

function continuarMedirTodo()
    DatosNombreArchivo = datosNombreArchivo(N, dimension, tipo, phi, crecimiento)
    
end

mysum(f, itr) = sum(f.(itr))
                                            
function dtvecinosDypo(Caja::caja{D}, Celdas::celdas{D}, Eventos::eventos, phi::Float64, crecimiento::Float64; Dt::Float64 = 0.5, modNombre::String="") where D
    DatosNombreArchivo = datosNombreArchivo(Caja.N, Caja.D, Caja.Ntipos, phi, crecimiento)
    NombreArchivoPO6T = "PO6T"*modNombre*DatosNombreArchivo*".csv"
            
    Ca, Ce, Ev = (deepcopy(Caja), deepcopy(Celdas), deepcopy(Eventos))
                                                                                
    tiempos = Float64[]
    POS6 = Array{Float64,1}[]
    POS12 = Array{Float64,1}[]
    AvgVNcols = Array{Float64,1}[]
    AvgCols = Array{Float64,1}[]
                               
    tiempoInicial = Ca.tiempoG    
    tiempo = 0.0
    Ncols = [1, 2, 3, 5, 10, 20, 50, 100]
    dt = 0.001
    while Ca.tiempoG < tiempoInicial + Dt
        
        Ca, Ce, Ev = evolucionardt(Ca, Ce, Ev, dt)
        Ca = mover_esferas(Ca)

        AvgVNcol = Float64[]
        AvgCol = Float64[]
        PO6 = Float64[]
        PO12 = Float64[]
        for n in Ncols
            
            # #sum([], init=0) julia 1.6
            # push!(AvgVNcol, sum(count.(x->(x.ncols >= n), Caja.vecinos)) / Caja.N)
            # push!(AvgCol, sum(sum.(x->if x.ncols >= n x.ncols else 0 end, Caja.vecinos)) / Caja.N)

            push!(AvgVNcol, mysum(x->x, count.(x->(x.ncols >= n), Ca.vecinos)) / Ca.N)
            push!(AvgCol, mysum(x->x, mysum.(x->if x.ncols >= n x.ncols else 0 end, Ca.vecinos)) / Ca.N)

            if D == 2
                PO = psi(Ca, 6, ncols=n), psi(Ca, 12, ncols=n)
                push!(PO6, PO[1])
                push!(PO12, PO[2])
            elseif D == 3
                push!(PO6, Q6(Ca, ncols=n))
                push!(PO12, 0.0)
            end
        end
        push!(AvgVNcols, AvgVNcol)
        push!(AvgCols, AvgCol)
        push!(POS6, PO6)
        push!(POS12, PO12)
        
        tiempo += dt
        tiempo = round(tiempo, digits=3)
        push!(tiempos, tiempo)

        if tiempo ∈ [0.002, 0.005, 0.01, 0.03, 0.05, 0.1, 0.3, 0.5]                 
            mkpath("GraficosVecinosdt")
            medirDistribucionVecinosDinamicos(Ca, Ce, Ev, phi, crecimiento, modNombre=modNombre*"dt_$tiempo"*"_")
            
            if D == 2
#                 graficaVecinos(Ca, Ce, tiempo)
#                 savefig("GraficosVecinosdt/GV"*modNombre*"dt_"*"$tiempo"*"_"*DatosNombreArchivo*".png")
                # savefig("GV"*modNombre*"$tiempo"*DatosNombreArchivo*".pdf")

                # graficamucz(Caja, Celdas)
                # savefig("GM"*modNombre*"$tiempo"*DatosNombreArchivo*".png")
                # savefig("GM"*modNombre*"$tiempo"*DatosNombreArchivo*".pdf")
            end
        end
    end
    Ncolss = fill(Ncols, length(tiempos))
    CSV.write(NombreArchivoPO6T, DataFrame(tiempo = tiempos, PO6 = POS6, PO12 = POS12, AvgC = AvgCols, AvgV = AvgVNcols, Ncols = Ncolss ))
end

function dtEvecinosDypo(Caja::caja{D}, Celdas::celdas{D}, Eventos::eventos, SqrtE::Float64, phi::Float64, crecimiento::Float64; Nt::Int64 = 500, modNombre::String="") where D
    DatosNombreArchivo = datosNombreArchivo(Caja.N, Caja.D, Caja.Ntipos, phi, crecimiento)
    NombreArchivoPO6T = "PO6TE"*modNombre*DatosNombreArchivo*".csv"
    
    Ca, Ce, Ev = (deepcopy(Caja), deepcopy(Celdas), deepcopy(Eventos))
                                                                                        
    VPO = VecinosParametrodeOrden(Ca)
    POotros = 0.0
    if D == 2
        POotros = psi6(Ca, VPO) 
    elseif D == 3
        POotros = Q6(Ca, VPO) 
    end
    CSV.write( "PO-Otros-"*modNombre*DatosNombreArchivo*".csv", DataFrame(POotros = POotros))
    tiempos = Float64[]
    nts = Int64[]
    POS6 = Array{Float64,1}[]
    POS12 = Array{Float64,1}[]
    AvgVNcols = Array{Float64,1}[]
    AvgCols = Array{Float64,1}[]
                
    tiempoInicial = Ca.tiempoG
    Ncols = [1, 2, 3, 5, 10, 20, 50, 100]
    dt = 1 / (SqrtE)
    nt = 0
        CSV.write("dtE"*DatosNombreArchivo*".csv", DataFrame(dt = dt, SqrtE = SqrtE))
    while nt < Nt
        
        Ca, Ce, Ev = evolucionardt(Ca, Ce, Ev, dt)
        Ca = mover_esferas(Ca)

        AvgVNcol = Float64[]
        AvgCol = Float64[]
        PO6 = Float64[]
        PO12 = Float64[]
        for n in Ncols
            
            # #sum([], init=0) julia 1.6
            # push!(AvgVNcol, sum(count.(x->(x.ncols >= n), Caja.vecinos)) / Caja.N)
            # push!(AvgCol, sum(sum.(x->if x.ncols >= n x.ncols else 0 end, Caja.vecinos)) / Caja.N)

            push!(AvgVNcol, mysum(x->x, count.(x->(x.ncols >= n), Ca.vecinos)) / Ca.N)
            push!(AvgCol, mysum(x->x, mysum.(x->if x.ncols >= n x.ncols else 0 end, Ca.vecinos)) / Ca.N)
        
                                                                              
            if D == 2
                PO = psi(Ca, 6, ncols=n), psi(Ca, 12, ncols=n)
                push!(PO6, PO[1])
                push!(PO12, PO[2])
            elseif D == 3
                push!(PO6, Q6(Ca, ncols=n))
                push!(PO12, 0.0)
            end
        end
        push!(AvgVNcols, AvgVNcol)
        push!(AvgCols, AvgCol)
        push!(POS6, PO6)
        push!(POS12, PO12)
        
        nt += 1        
        push!(tiempos, Ca.tiempoG - tiempoInicial)
        push!(nts, nt)

        if nt ∈ [500, 1000]
            medirDistribucionVecinosDinamicos(Ca, Ce, Ev, phi, crecimiento, modNombre=modNombre*"nt_$(nt)"*"_")
            
            if D == 2
#                 mkpath("GraficosVecinosdtE") 
#                 graficaVecinos(Ca, Ce, )
#                 savefig("GraficosVecinosdtE/GV"*modNombre*"_nt_$(nt)"*"_"*DatosNombreArchivo*".png")
                # savefig("GV"*modNombre*"$tiempo"*DatosNombreArchivo*".pdf")

                # graficamucz(Caja, Celdas)
                # savefig("GM"*modNombre*"$tiempo"*DatosNombreArchivo*".png")
                # savefig("GM"*modNombre*"$tiempo"*DatosNombreArchivo*".pdf")
            end
        end
    end
    Ncolss = fill(Ncols, length(tiempos))
    CSV.write(NombreArchivoPO6T, DataFrame(tiempo = tiempos, nts = nts, PO6 = POS6, PO12 = POS12, AvgC = AvgCols, AvgV = AvgVNcols, Ncols = Ncolss ))
end

function CpVecinosDyPO(Caja::caja{D}, Celdas::celdas{D}, Eventos::eventos, phi::Float64, crecimiento::Float64; modNombre::String="") where D
    DatosNombreArchivo = datosNombreArchivo(Caja.N, Caja.D, Caja.Ntipos, phi, crecimiento)
    NombreArchivoPO6T = "PO6T"*modNombre*DatosNombreArchivo*".csv"
   
    tiempos = Float64[]
    POS = Float64[]
    cols = Float64[]

    tiempoInicial = Caja.tiempoG
    factor = Int64[10, 20, 50, 100, Caja.N/5, Caja.N/2, Caja.N]
    cont = 1
    while sum(Caja.numeroColisiones) <= Caja.N^2 + 2
        Caja, Celdas, Eventos = proximo_evento(Caja, Celdas, Eventos)

        Caja = mover_esferas(Caja)
        if D == 2
            push!(POS, psi6(Caja))
        elseif D == 3
            push!(POS, Q6(Caja))
        end
        push!(tiempos, Caja.tiempoG - tiempoInicial)
        push!(cols, round(sum(Caja.numeroColisiones)/Caja.N, digits=3))

        if factor[cont]*Caja.N -2 <= sum(Caja.numeroColisiones) <= factor[cont]*Caja.N + 2
            FC = factor[cont]
            cont += 1
            medirDistribucionVecinosDinamicos(Caja, Celdas, Eventos, phi, crecimiento, modNombre=modNombre*"FC$FC")
            
            if D == 2
                graficaVecinos(Caja, Celdas)
                savefig("GV"*modNombre*"$FC"*DatosNombreArchivo*".png")
                savefig("GV"*modNombre*"$FC"*DatosNombreArchivo*".pdf")

                graficamucz(Caja, Celdas)
                savefig("GM"*modNombre*"$FC"*DatosNombreArchivo*".png")
                savefig("GM"*modNombre*"$FC"*DatosNombreArchivo*".pdf")
            end
            if cont == 7
                break
            end
        end
    end
    CSV.write(NombreArchivoPO6T, DataFrame(PO = POS, tiempo = tiempos, AvgC = cols ))
end

function medirSoloDifusion(N::Int64, tipo::Int64, D::Int64, crecimiento::Float64, phi::Float64; propGrCh=1.4)
       # tdif en args
    println("Creciendo particulas con phi = $phi , con $tipo radios en $D D")
    Caja, Celdas, Eventos = CrecerDiscosAPhi(N, tipo, D,  crecimiento, phi, propGC=propGrCh)
#     Eventos = eventos2(Caja.tiempoG, Caja, Celdas)
    println("Termalizando")
    Caja, Celdas, Eventos = evolucionarFlu(Caja, Celdas, Eventos, 0.1)
    Caja = reiniciar(Caja)
    Caja = mover_esferas(Caja)

    Caja, Celdas, Eventos = evolucionarFlu(Caja, Celdas, Eventos, 9.0, addvecinos=false)
    Caja = reiniciar(Caja)
                                                                                                
    if D == 2
        tdif = 50.0
    elseif D == 3
        tdif = 50.0
    end
                                                                                                
    println("Midiendo difusion")
    medirDifusion(Caja, Celdas, Eventos, crecimiento, phi, tdif)
end

function medirTodo(N::Int64, tipo::Int64, D::Int64, crecimiento::Float64, phi::Float64; propGrCh=1.4)
    DatosNombreArchivo = datosNombreArchivo(N, D, tipo, phi, crecimiento)
#     dirs = ["Crecido", "Est", "GrF"]
                                                              
    println("Creciendo particulas con phi = $phi , con $tipo radios en $D D")
    Caja, Celdas, Eventos = CrecerDiscosAPhi(N, tipo, D,  crecimiento, phi, propGC=propGrCh)
#     Eventos = eventos2(Caja.tiempoG, Caja, Celdas)
    println("Termalizando")
    Caja, Celdas, Eventos = evolucionarFlu(Caja, Celdas, Eventos, 0.1)
    Caja = mover_esferas(Caja)
    # saveAll(Caja, Eventos, dirs[1])

#     graficaST(Caja, Celdas)
#     savefig(DatosNombreArchivo*".png")
#     D == 2 ? savefig(DatosNombreArchivo*".pdf") : nothing
    Caja = reiniciar(Caja)
    
    saveXV(Caja, "I")
    SqrtEi = medirEnergia(Caja, Celdas, Eventos, phi, crecimiento, modNombre="I")
    dtEvecinosDypo(Caja, Celdas, Eventos, SqrtEi, phi, crecimiento, Nt=1000, modNombre="I")
    
#     println("Obteniendo g(r)")
#     if tipo == 1 
#         medirDistribucionRadialP(Caja, Celdas, Eventos, phi, crecimiento, modNombre="I")
#     elseif tipo == 2
#         medirDistribucionRadialBinariop(Caja, Celdas, Eventos, phi, crecimiento, modNombre="I")
#     end

#      if rapido 
                                                                                                        
#                                                                                                         end
                                                                                                    
    Caja, Celdas, Eventos = evolucionarFlu(Caja, Celdas, Eventos, 9.0, addvecinos=false)
    Caja = mover_esferas(Caja)
    Caja = reiniciar(Caja)

    # println("Midiendo difusion")
    # medirDifusion(Caja, Celdas, Eventos, phi, crecimiento, 1.0)
    # Caja = reiniciar(Caja)
    Caja = mover_esferas(Caja)

    saveXV(Caja, "F")
    SqrtEf = medirEnergia(Caja, Celdas, Eventos, phi, crecimiento, modNombre="F")
    dtEvecinosDypo(Caja, Celdas, Eventos, SqrtEf, phi, crecimiento, Nt=1000, modNombre="F")

#     println("Obteniendo parametro de orden")
#     medirParametroOrdenOrientacional(Caja, Celdas, Eventos, phi, crecimiento, modNombre="F")
#     println("Obteniendo g(r)")
#     if tipo == 1 
#         medirDistribucionRadialP(Caja, Celdas, Eventos, phi, crecimiento, modNombre="F")
#     elseif tipo == 2
#         medirDistribucionRadialBinarioP(Caja, Celdas, Eventos, phi, crecimiento, modNombre="F")
#     end
    # saveAll(Caja, Eventos, dirs[3])

    # println("Obteniendo lineas de tiempo")
    # medirLineasDeTiempo(Caja, Celdas, Eventos, phi, crecimiento, 1.0)
    # println("Termino")
end
                                                                    
function medirSoloSVD(N::Int64, tipo::Int64, D::Int64, crecimiento::Float64, phi::Float64; propGrCh=1.4)
    DatosNombreArchivo = datosNombreArchivo(N, D, tipo, phi, crecimiento)
#     dirs = ["Crecido", "Est", "GrF"]
                                                                
    println("Creciendo particulas con phi = $phi , con $tipo radios en $D D")
    Caja, Celdas, Eventos = CrecerDiscosAPhi(N, tipo, D,  crecimiento, phi, propGC=propGrCh)
    Eventos = eventos(Caja.tiempoG, Caja, Celdas)
    println("Termalizando")
    Caja, Celdas, Eventos = evolucionarFlu(Caja, Celdas, Eventos, 0.2)
    Caja = mover_esferas(Caja)
    SqrtEi = medirEnergia(Caja, Celdas, Eventos, phi, crecimiento, modNombre="I")
    # saveAll(Caja, Eventos, dirs[1])

    graficaST(Caja, Celdas)
    savefig(DatosNombreArchivo*".png")
    D == 2 ? savefig(DatosNombreArchivo*".pdf") : nothing
    Caja = reiniciar(Caja)

    dtvecinosDypo(Caja, Celdas, Eventos, phi, crecimiento, Dt=0.5, modNombre="I")
    dtEvecinosDypo(Caja, Celdas, Eventos, SqrtEi, phi, crecimiento, Nt=1000, modNombre="I")
       
    Caja, Celdas, Eventos = evolucionarFlu(Caja, Celdas, Eventos, 9.0, addvecinos=false)
    Caja = mover_esferas(Caja)
    Caja = reiniciar(Caja)

    SqrtEf = medirEnergia(Caja, Celdas, Eventos, phi, crecimiento, modNombre="F")
    dtvecinosDypo(Caja, Celdas, Eventos, phi, crecimiento, Dt=0.5, modNombre="F")
    dtEvecinosDypo(Caja, Celdas, Eventos, SqrtEf, phi, crecimiento, Nt=1000, modNombre="F")
end

function medirParametroOrdenOrientacionalX(Caja::caja{D}, Celdas::celdas{D}, Eventos::eventos, phi::Float64, crecimiento::Float64, x=400.0; modNombre::String="") where D
    DatosNombreArchivo = datosNombreArchivo(Caja.N, Caja.D, Caja.Ntipos, phi, crecimiento)
    NombreArchivoPO = "PO"*modNombre*DatosNombreArchivo*".csv"

    if D == 2
        NombreArchivoPOX = "POX"*modNombre*DatosNombreArchivo*".csv"
        O = 0:0.01:x
        POX = zeros(length(O))
        for i in 1:length(O)
            POX[i] = psi(Caja, O[i])
        end
        CSV.write(NombreArchivoPOX, DataFrame(orden = O, psix = POX))
        
        plot(title = "Con $Caja.ntipos radios, $Caja.N discos",
            xlabel = "Orden",  ylabel = "psi x", label = "phi = $phi")
        plot!(O, POX, label = "")
        savefig("pox-"*modNombre*DatosNombreArchivo*".png")
    end
end                                                                                            

function medirSimetrias(N::Int64, tipo::Int64, D::Int64, crecimiento::Float64, phi::Float64; propGCrn=[1.4, 2.0])
    DatosNombreArchivo = datosNombreArchivo(N, D, tipo, phi, crecimiento)
                                                              
    println("Creciendo particulas con phi = $phi , con $tipo radios en $D D")
    Caja, Celdas, Eventos = CrecerDiscosAPhi(N, tipo, D,  crecimiento, phi, propGCrn=propGCrn)
    Eventos = eventos(Caja.tiempoG, Caja, Celdas)
    println("Termalizando")
    Caja, Celdas, Eventos = evolucionarFlu(Caja, Celdas, Eventos, 0.1)
    Caja = mover_esferas(Caja)
    Caja = reiniciar(Caja)
    
    saveXV(Caja, "I")
    SqrtEi = medirEnergia(Caja, Celdas, Eventos, phi, crecimiento, modNombre="I")
#     dtEvecinosDypo(Caja, Celdas, Eventos, SqrtEi, phi, crecimiento, Nt=1000, modNombre="I")
#     medirParametroOrdenOrientacionalX(Caja, Celdas, Eventos, phi, crecimiento, 500.0,  modNombre="I")
                                                                                                    
    Caja, Celdas, Eventos = evolucionarFlu(Caja, Celdas, Eventos, 9.0, addvecinos=false)
    Caja = mover_esferas(Caja)
    Caja = reiniciar(Caja)

    saveXV(Caja, "F")
    SqrtEf = medirEnergia(Caja, Celdas, Eventos, phi, crecimiento, modNombre="F")
#     dtEvecinosDypo(Caja, Celdas, Eventos, SqrtEf, phi, crecimiento, Nt=1000, modNombre="F")
#     medirParametroOrdenOrientacionalX(Caja, Celdas, Eventos, phi, crecimiento, 500.0,  modNombre="F")

    graficaST(Caja, Celdas)
    savefig(DatosNombreArchivo*".png")
    savefig(DatosNombreArchivo*".pdf")
end
                                                                                                
function medirCSimetrias(N::Int64, tipo::Int64, D::Int64, crecimiento::Float64, n::Int64; propGCrn=[1.4, 2.0])
    println("Creciendo particulas con n = $n pasos, con $tipo radios en $D D")
                                                                                                    
    Caja, Celdas, Eventos = crecer(N, tipo, D, crecimiento, n; propGCrn=propGCrn, v=5.0)
    Caja = parar_crecimiento(Caja)
    Caja = mover_esferas(Caja)
    @show phi = packingFraction(Caja)
    DatosNombreArchivo = datosNombreArchivo(N, D, tipo, phi, crecimiento)
                                                              
    Eventos = eventos(Caja.tiempoG, Caja, Celdas)
    println("Termalizando")
    Caja, Celdas, Eventos = evolucionarFlu(Caja, Celdas, Eventos, 0.1)
    Caja = mover_esferas(Caja)
    Caja = reiniciar(Caja)
    
    saveXV(Caja, "I")
    SqrtEi = medirEnergia(Caja, Celdas, Eventos, phi, crecimiento, modNombre="I")
#     dtEvecinosDypo(Caja, Celdas, Eventos, SqrtEi, phi, crecimiento, Nt=1000, modNombre="I")
#     medirParametroOrdenOrientacionalX(Caja, Celdas, Eventos, phi, crecimiento, 500.0,  modNombre="I")
                                                                                                    
    Caja, Celdas, Eventos = evolucionarFlu(Caja, Celdas, Eventos, 9.0, addvecinos=false)
    Caja = mover_esferas(Caja)
    Caja = reiniciar(Caja)

    saveXV(Caja, "F")
    SqrtEf = medirEnergia(Caja, Celdas, Eventos, phi, crecimiento, modNombre="F")
#     dtEvecinosDypo(Caja, Celdas, Eventos, SqrtEf, phi, crecimiento, Nt=1000, modNombre="F")
#     medirParametroOrdenOrientacionalX(Caja, Celdas, Eventos, phi, crecimiento, 500.0,  modNombre="F")

    graficaSP(Caja, Celdas)
    savefig(DatosNombreArchivo*".png")
    savefig(DatosNombreArchivo*".pdf")
end
        
function proximo_eventoCn1(Caja::caja{D}, Celdas::celdas{D}, Eventos::eventos; vecinos::Bool = false, dtV = 0.0) where D
    tiempoMin = minimum(Eventos.tiemposCyT)
    n1, j = findfirst(Eventos.tiemposCyT, tiempoMin)
    
    Caja.tiempoG = tiempoMin
    n2 = Eventos.eventoCyT[n1, j]
    coln1 = 0
    if j == 1 #colision
        if n2 > 0 #con esfera
            Caja = mover_esfera(n1, Caja)
            Caja = mover_esfera(n2, Caja)
            Caja = colision(n1, n2, Caja)
            if vecinos
                Caja = agregar_vecinos(n1, n2, Caja)
                if dtV != 0.0
                    Caja = quitar_vecinos!(Caja, dtV)
                end
            end
            coln1 = n1
            
            Eventos = tiempo_transferencia_nuevo(n1, Caja, Celdas, Eventos)
            Eventos = tiempo_transferencia_nuevo(n2, Caja, Celdas, Eventos)
            
            Eventos = tiempos_colision_nuevos((n1,n2), Caja, Celdas, Eventos)
            Eventos = recalcular_tiempo_colision((n1,n2), Caja, Celdas, Eventos)
        else #paredes
            Caja = colision_pared(n1, n2, Caja)
            
            Eventos = tiempo_transferencia_nuevo(n1, Caja, Celdas, Eventos)
            Eventos = tiempos_colision_nuevos(n1, Caja, Celdas, Eventos)
        end
    else #tranferencia
        Celdaold = Celdas.presenteEnCelda[n1]
        Caja, Celdas = transferenciaCelda(n1, n2, Caja, Celdas)
        
        Eventos = tiempo_transferencia_nuevo(n1, Caja, Celdas, Eventos)
        Eventos = tiempo_colision_PTransferencia(n1, Celdaold, Caja, Celdas, Eventos)
    end
    Caja, Celdas, Eventos, coln1
end

function dtAverageMeassure(Caja, Celdas, Eventos)# Después de crecer con número de pasos
    contT = 0
    t0 = Caja.tiempoG
    dtA = 0.0
    tsC = zeros(Caja.N)
    criterioAlcanzado = false
    while !criterioAlcanzado
        
        Caja, Celdas, Eventos, n1 = proximo_eventoCn1(Caja, Celdas, Eventos)
        if n1 != 0
            if tsC[n1] == 0.0
        
                tsC[n1] = Caja.tiempoG
        
            elseif tsC[n1] != 0.0
                dtA += Caja.tiempoG - tsC[n1]
                contT += 1
                tsC[n1] = Caja.tiempoG
            end
        end
        if contT > 2N
            return dtA / contT
        end
    end
end
         
function dtAtasque(Caja, Celdas, Eventos; dt::Float64=1e-10)
    contT = 0
    t0 = Caja.tiempoG
    dtA = 0.0
    tsC = zeros(Caja.N)
    criterioAlcanzado = false
    while !criterioAlcanzado
        
        Caja, Celdas, Eventos, n1 = proximo_eventoCn1(Caja, Celdas, Eventos)
        if n1 != 0
            if tsC[n1] == 0.0
        
                tsC[n1] = Caja.tiempoG
        
            elseif tsC[n1] != 0.0
                dtA += Caja.tiempoG - tsC[n1]
                contT += 1
                tsC[n1] = Caja.tiempoG
            end
        end
        if contT > 2*Caja.N
            if dt > dtA / contT 
                @show dtA / contT
                criterioAlcanzado = true
            else
                contT = 0
                dtA = 0.0
            end
        end
    end
    Caja, Celdas, Eventos
end                                                                                       
        
function crecerAJam(N::Int64, tipo::Int64, D::Int64, crecimiento::Float64; dtTarget::Float64 = 5e-10, propGCrn::Array{Float64,1}=[1.4,2.0], v::Float64=5.0)
    Caja, Celdas, Eventos = condiciones_iniciales(N, tipo, D, crecimiento, prop = propGCrn, vm = v)
    Caja, Celdas, Eventos = dtAtasque(Caja, Celdas, Eventos, dt=dtTarget)
    
    Caja = parar_crecimiento(Caja)
    Caja = mover_esferas(Caja)             
    Eventos = eventos(Caja.tiempoG, Caja, Celdas)
    phi = packingFraction(Caja)
    Caja, Celdas, Eventos, phi
end                                                                                        
                                                                                                
function medirSimetriasAtascadas(N::Int64, tipo::Int64, D::Int64, crecimiento::Float64, dtTarget::Float64; propGCrn=[1.4, 2.0])
    println("Creciendo particulas hasta alcanzar el promedio entre colisiones de dt = $dt, con $tipo radios con proporciones [r,p] = $propGCrn en $D D")
                                                                                                    
    Caja, Celdas, Eventos, phi = crecerAJam(N, tipo, D, crecimiento, dtTarget=dt, propGCrn=propGCrn)
    @show phi  
                                                                                                                        
    DatosNombreArchivo = datosNombreArchivo(N, D, tipo, phi, crecimiento)
                           
    saveXV(Caja, "I")
    SqrtEi = medirEnergia(Caja, Celdas, Eventos, phi, crecimiento, modNombre="I")
                                                                                                    

    graficaSP(Caja, Celdas)
    savefig(DatosNombreArchivo*".pdf")
end
                                                                                    


# Caja, Celdas, Eventos = condiciones_inicialesTrian(propRadioMax=0.995)
# # phimaax = 0.899

# Caja, Celdas, Eventos = reconstruir(nombreDir::String)
# ["Crecido", "Est", "GrF"]


# medirCasiTodo(10_000, 1, 2, 8e-3, 0.82)

# medirCasiTodo(1_000, 2, 2, 8e-3, 0.82, propGrCH = 1.4)

# medirCasiTodo(1_500, 1, 3, 8e-3, 0.52)