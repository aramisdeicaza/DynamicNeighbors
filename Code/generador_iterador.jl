function mover_esfera(n1::Int64, Caja::caja{D}) where D
    s1 = Caja.esferas[n1]
    Δt = Caja.tiempoG - s1.t
    Caja.esferas[n1].x = s1.x + s1.v*Δt
    if s1.c != 0.0
        Caja.esferas[n1].r = s1.r + s1.c*(Caja.tiempoG - s1.t)
        # Caja.esferas[n1].r = s1.c*Caja.tiempoG
    end 
    Caja.esferas[n1].t = Caja.tiempoG
    Caja
end

function mover_esferas(Caja::caja{D}) where D
    for n1 in 1:Caja.N
        Caja = mover_esfera(n1, Caja)
    end
    Caja
end


function mover_esferas(Caja::caja{D}, Celdas::celdas, Eventos::eventos) where D
    for n1 in 1:Caja.N
        Caja = mover_esfera(n1, Caja)
    end
    Eventos = eventos(Caja.tiempoG, Caja, Celdas)
    Caja, Celdas, Eventos
end

function colision(n1::Int64, n2::Int64, Caja::caja{D}) where D
    s1 = Caja.esferas[n1]
    s2 = Caja.esferas[n2]
    
    v1 = s1.v
    v2 = s2.v

    Δp = s1.x - s2.x
    Δp = distanciaPorComponente(Δp, Caja)
    
    m1 = s1.m
    m2 = s2.m

    n = normalize(Δp)  
    v1n = (v1 ⋅ n) * n  
    v2n = (v2 ⋅ n) * n

    propMasa = (1 / ((m1 / m2) + 1))
    
    if s1.c != 0
        vc = n*(s1.c + s2.c)
    else
        vc = zeros(D)
    end
    v1f = v1 + 2 * propMasa * (v2n - v1n + vc)
    v2f = v2 + 2 * (1 - propMasa) * (v1n - v2n - vc)

    Caja.esferas[n1].v = v1f
    Caja.esferas[n2].v = v2f

    Caja
end

function colision_pared(n1::Int64, colconpared::Int64, Caja::caja{D}) where D
        #Falta agregar empuje del crecimiento
    s1 = Caja.esferas[n1]
    
    i = -floor(Int, colconpared / 2)
    r = SMatrix{D,D}(ones(D,D) - 2I)
    
    Caja.esferas[1].v = s1.v.*r[i,:]
    
    Caja
end

function transferenciaCelda(n1::Int64, transferencia::Int64, Caja::caja{D}, Celdas::celdas{D}) where D
    ijk = copy(Celdas.presenteEnCelda[n1])
    
    for l in 1:length(Celdas.residentes[ijk])
        if Celdas.residentes[ijk][l] == n1
            deleteat!(Celdas.residentes[ijk], l)
#             error("No se encuentra $n1 en celda $ijk")
            break
        end
    end
    
    correcionPosiciones = zeros(D)
    cont = 0
    for i in 1:D
        for j in 1:2
            cont += 1
            if j == 1 && transferencia  == cont
                ijk[i] += 1
                if ijk[i] == Celdas.m + 1
                    ijk[i] = 1
                    correcionPosiciones[i] = -Caja.L
                end
                break
            elseif j == 2 && transferencia == cont
                ijk[i] -= 1
                if ijk[i] == 0
                    ijk[i] = Celdas.m 
                    correcionPosiciones[i] = Caja.L
                end
                break
            end
        end
    end
       
    Caja.esferas[n1].x += correcionPosiciones
    Caja.transferenciasPeriodicas[n1] += -correcionPosiciones
    
    push!(Celdas.residentes[ijk], n1)
    Celdas.presenteEnCelda[n1] = ijk
     
    Caja, Celdas
end

function agregar_vecinos(n1::Int64, n2::Int64, Caja::caja{D}) where D
    vecinon1 = 0
    for i in 1:Caja.numeroVecinos[n1]
        if Caja.vecinos[n1][i].vecino == n2
            vecinon1 = i
            break
        end
    end
    if vecinon1 == 0
        insert!(Caja.vecinos[n1], vecinosDinamicos(n2, Caja.tiempoG, 1, Caja.esferas[n2].tipo))
        insert!(Caja.vecinos[n2], vecinosDinamicos(n1, Caja.tiempoG, 1, Caja.esferas[n1].tipo))

        Caja.numeroVecinos[n1] += 1
        Caja.numeroVecinos[n2] += 1
        if Caja.Ntipos == 2
            Caja.tipoDeVecinos[n1, Caja.esferas[n2].tipo] += 1 
            Caja.tipoDeVecinos[n2, Caja.esferas[n1].tipo] += 1
        end
    else
        vecinon2 = 0
        for i in 1:Caja.numeroVecinos[n2]
            if Caja.vecinos[n2][i].vecino == n1
                vecinon2 = i
                break
            end
        end
        Caja.vecinos[n1][vecinon1].tiempo = Caja.tiempoG
        Caja.vecinos[n1][vecinon1].ncols += 1
    
        Caja.vecinos[n2][vecinon2].tiempo = Caja.tiempoG
        Caja.vecinos[n2][vecinon2].ncols += 1
    end
    Caja.numeroColisiones[n1] += 1
    Caja.numeroColisiones[n2] += 1
    Caja
end

function quitar_vecinos!(Caja::caja{D}, dtV::Float64) where D
    for n1 in 1:Caja.N
        sort!(Caja.vecinos[n1])
        n = count(x->(Caja.tiempoG - x.tiempo > dtV), Caja.vecinos[n1])
        if n != 0 ncolsn = sum(x->x.ncols, Caja.vecinos[n1][1:n]) else ncolsn = 0 end
        Caja.vecinos[n1] = Caja.vecinos[n1][n+1:end]
        Caja.numeroVecinos[n1] -= n
        Caja.numeroColisiones[n1] -= ncolsn
    end
    Caja
end

function condiciones_iniciales(N::Int64, ntipos::Int64, D::Int64, Crecimiento::Float64; prop::Array{Float64,1} = [1.4, 2.0], vm::Float64 = 5.0)
    Caja = caja(N, ntipos, D, Crecimiento, propGCrn = prop, vmax = vm)
    Celdas = celdas(Caja)
    Eventos = eventos(Caja, Celdas)
    Caja, Celdas, Eventos
end

# function condiciones_inicialesTrian(propRadioMax::Float64=0.995)
#     N, r, pos, phiId = triangular(1, 184, 186)
#     r = propRadioMax*r
#     phi = N*pi*r^2
#     @show phi, phiId
#     Caja = caja(N, r, pos)
#     Celdas = celdas(Caja)
#     Eventos = eventos(Caja, Celdas)
#     Caja, Celdas, Eventos
# end

# function reconstruir(nombreDir::String)
#     Caja = caja(nombreDir)
#     Celdas = celdas(Caja)
#     Eventos = eventos(nombreDir)
#     Caja, Celdas, Eventos
# end


import Base: findfirst

function findfirst(X::Array{Float64,2}, a::Float64)# where T<:Real
    for i in 1:length(X[:,1])
        for j in 1:2
            if X[i, j] == a
                return i, j
            end
        end
    end
    return 0
end

function proximo_evento(Caja::caja{D}, Celdas::celdas{D}, Eventos::eventos; agregarvecinos::Bool = true, dtV = 0.0) where D
    tiempoMin = minimum(Eventos.tiemposCyT)
    n1, j = findfirst(Eventos.tiemposCyT, tiempoMin)
    
    Caja.tiempoG = tiempoMin
    n2 = Eventos.eventoCyT[n1, j]

    if j == 1 #colision
        if n2 > 0 #con esfera
            Caja = mover_esfera(n1, Caja)
            Caja = mover_esfera(n2, Caja)
            Caja = colision(n1, n2, Caja)
            if agregarvecinos
                Caja = agregar_vecinos(n1, n2, Caja)
                if dtV != 0.0
                    Caja = quitar_vecinos!(Caja, dtV)
                end
                # checar con dtv
            end
            
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
    Caja, Celdas, Eventos
end

function crecer(N::Int64, tipo::Int64, D::Int64, crecimiento::Float64, n::Int64; propGCrn::Array{Float64,1}=[1.4,2.0], v::Float64=5.0)
    Caja, Celdas, Eventos = condiciones_iniciales(N, tipo, D, crecimiento, prop = propGCrn, vm = v)
    prints = 1
    for i in 1:n
        if prints/5 < i/n && prints <= 5
            println("Avance de crecimiento a $prints / 5")
            prints += 1
        end
        Caja, Celdas, Eventos = proximo_evento(Caja, Celdas, Eventos, agregarvecinos = false)
    end
    Caja, Celdas, Eventos
end

function reiniciar(Caja::caja{D}) where D
    for n1 in 1:Caja.N
        Caja.numeroVecinos[n1] = 0
        Caja.numeroColisiones[n1] = 0
        Caja.vecinos[n1] = []
        Caja.transferenciasPeriodicas[n1] = zeros(D)
    end
    Caja
end

function parar_crecimiento(Caja::caja{D}) where D
    Caja = reiniciar(Caja)
    for n1 in 1:Caja.N
        s1 = Caja.esferas[n1]
        Caja.esferas[n1].r = s1.r + (Caja.tiempoG - s1.t)*s1.c
        Caja.esferas[n1].c = 0.0
    end
    Caja
end

function CrecerDiscosAPhi(N::Int64, tipo::Int64, D::Int64, crecimiento::Float64, ϕ::Float64; propGCrn::Array{Float64,1} = [1.4, 2.0], v::Float64=5.0)
    Caja, Celdas, Eventos = condiciones_iniciales(N, tipo, D, crecimiento, prop = propGCrn, vm = v)
    tPhi = 0
    if D == 2
        if tipo == 1
            tPhi = (ϕ / (N*pi*crecimiento^2))^(1/2)
        elseif tipo == 2
            N1, N2 = contarGrCh(Caja.esferas)
            tPhi = (ϕ / (pi*N1*crecimiento^2 + pi*N2*(Caja.propGCrn[1]*crecimiento)^2))^(1/2)
        end
    elseif D == 3
        tPhi = ( (3/4)*ϕ / (N*pi*crecimiento^3))^(1/3)
    end
    
    if Celdas.l < 2*tPhi*Caja.propGCrn[1]*crecimiento
        @show Celdas.l, 2*tPhi*Caja.propGCrn[1]*crecimiento
        error("radio sera grande")
    end
            
    prints = 1
    cont = 0
    while Caja.tiempoG < tPhi
        if prints*tPhi/5 < Caja.tiempoG
            println("Avance de crecimiento a $prints / 5")
            prints += 1
        end
        Caja, Celdas, Eventos = proximo_evento(Caja, Celdas, Eventos, agregarvecinos=false)
    end
   
    Caja = parar_crecimiento(Caja)
    Caja = mover_esferas(Caja)
    Eventos = eventos(Caja.tiempoG, Caja, Celdas)
    Caja, Celdas, Eventos
end

function evolucionarNumeroPasos(Caja::caja{D}, Celdas::celdas{D}, Eventos::eventos, pasos::Int64) where D
    for i in 1:pasos
        Caja, Celdas, Eventos = proximo_evento(Caja, Celdas, Eventos)
    end
    Caja, Celdas, Eventos
end

function evolucionarFlu(Caja::caja{D}, Celdas::celdas{D}, Eventos::eventos, k::Float64; addvecinos::Bool=true) where D
    np::Int64 = ceil(k*Caja.N^2)
    for i in 1:np
        if k > 0.2 && i == floor(np/2)
            println("Avance de dinamica a la mitad")
        end
        Caja, Celdas, Eventos = proximo_evento(Caja, Celdas, Eventos, agregarvecinos=addvecinos)

    end
    Caja, Celdas, Eventos
end

function evolucionardt(Caja::caja{D}, Celdas::celdas{D}, Eventos::eventos, dt::Float64; addvecinos::Bool=true) where D
    t = Caja.tiempoG 
    while Caja.tiempoG < t + dt
        Caja, Celdas, Eventos = proximo_evento(Caja, Celdas, Eventos, agregarvecinos=addvecinos)
    end
    Caja, Celdas, Eventos
end

function date()
    Date = string(Dates.today())*"-"*string(Dates.hour(now()))*"h-"*string(Dates.minute(now()))*"m"
end

function datosNombreArchivo(N::Int64, dimension::Int64, tipo::Int64, phi::Float64, crecimiento::Float64)
    DatosNombreArchivo = "$dimension"*"D-"*"N_$N"*"-t_$tipo"*"-phi_$phi"*"-CR_$crecimiento"*"_"*date()
end

function dfxrt(Caja::caja{2})
    dfxrt = DataFrame(x1 = Caja.esferas[1].x[1], x2 = Caja.esferas[1].x[2], r = Caja.esferas[1].r, tipo = Caja.esferas[1].tipo)
    for i in 2:Caja.N
        append!(dfxrt, DataFrame( x1 = Caja.esferas[i].x[1], x2 = Caja.esferas[i].x[2], r = Caja.esferas[i].r, tipo = Caja.esferas[i].tipo))
    end
    dfxrt
end

function dfxrt(Caja::caja{3})
    dfxrt = DataFrame(x1 = Caja.esferas[1].x[1], x2 = Caja.esferas[1].x[2], x3 = Caja.esferas[1].x[3], 
                                            r = Caja.esferas[1].r, tipo = Caja.esferas[1].tipo)
    for i in 2:Caja.N
        append!(dfxrt, DataFrame(x1 = Caja.esferas[i].x[1], x2 = Caja.esferas[i].x[2], x3 = Caja.esferas[i].x[3], 
                                                    r = Caja.esferas[i].r, tipo = Caja.esferas[i].tipo))
    end
    dfxrt
end
                                    
function dfv(Caja::caja{2})
    dfv = DataFrame(v1 = Caja.esferas[1].v[1], v2 = Caja.esferas[1].v[2])
    for i in 2:Caja.N
        append!(dfv, DataFrame(v1 = Caja.esferas[i].v[1], v2 = Caja.esferas[i].v[2]))
    end
    dfv
end

function dfv(Caja::caja{3})
    dfv = DataFrame(v1 = Caja.esferas[1].v[1], v2 = Caja.esferas[1].v[2], v3 = Caja.esferas[1].v[3])
    for i in 2:Caja.N
        append!(dfv, DataFrame(v1 = Caja.esferas[i].v[1], v2 = Caja.esferas[i].v[2], v3 = Caja.esferas[i].v[3]))
    end
    dfv
end                                    

function saveXV(Caja::caja{D}, parte::String) where D
    #DataFrame no pone los Float64 bien
    phi = packingFraction(Caja)
    phi = round(phi, digits = 4)
    
    DatosNombreArchivo = datosNombreArchivo(Caja.N, Caja.D, Caja.Ntipos, phi, Caja.crecimientos[1])
    dir = parte*"-"*DatosNombreArchivo
    mkdir(dir)
       
#     NAEsferas = "Esferas.csv"
#     NAEventos = "Eventos.csv"
    NAXRT = "XRT.csv"
    NAV = "V.csv"
    NADatos = "Datos.csv"
#     NACaja = "Caja.csv"
    #fieldnames(typeof(Caja))
    
    if Caja.Ntipos == 1
        df = DataFrame(N = Caja.N, Ntipos = Caja.Ntipos, D = Caja.D, tiempoG = Caja.tiempoG, 
            propGCrn = Caja.propGCrn, cajon = Caja.caja, L = Caja.L, phi = phi, crecimiento = Caja.crecimientos[1])
    elseif Caja.Ntipos == 2
        N1, N2 = contarGrCh(Caja.esferas)
        df = DataFrame(N = Caja.N, Ntipos = Caja.Ntipos, N1 = N1, N2 = N2, D = Caja.D, tiempoG = Caja.tiempoG, 
            propGCrn = Caja.propGCrn, cajon = Caja.caja, L = Caja.L, phi = phi, crecimiento = Caja.crecimientos[1])   
    elseif Caja.Ntipos > 2
         error("Aun no implementado")                                   
    end
                                                
    CSV.write(dir*"/"*NADatos, df)

    CSV.write(dir*"/"*NAXRT, dfxrt(Caja))
    CSV.write(dir*"/"*NAV, dfv(Caja))

#     df1 = DataFrame(Caja.esferas)
#     CSV.write(dir*"/"*NAEsferas, df1)
        
    #quiza sea mejor reconstruirlo
#     df2 = DataFrame(tiemposC = Eventos.tiemposCyT[:,1], tiemposT = Eventos.tiemposCyT[:,2], 
#         eventosC = Eventos.eventoCyT[:,1], eventosT = Eventos.eventoCyT[:,2])
#     CSV.write(dir*"/"*NAEventos, df2)
        
    #agregar a struct de vecinosdinamicos?
#     df3 = DataFrame(transferenciaPeriodicas = Caja.transferenciasPeriodicas, numeroVecinos = Caja.numeroVecinos, 
#     numeroColisiones = Caja.numeroColisiones, vecinos = Caja.vecinos)
#     CSV.write(dir*"/"*NACaja, df3)
end
