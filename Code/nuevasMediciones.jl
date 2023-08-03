function medirnSoloDifusion(N::Int64, tipo::Int64, D::Int64, crecimiento::Float64, phi::Float64; propGrCh=1.4)
       # tdif en args
    println("Creciendo particulas con phi = $phi , con $tipo radios en $D D")
    Caja, Celdas, Eventos = CrecerDiscosAPhi(N, tipo, D,  crecimiento, phi, propGC=propGrCh)
    Eventos = eventos(Caja.tiempoG, Caja, Celdas)
    
    println("Termalizando")
    Caja, Celdas, Eventos = evolucionarFlu(Caja, Celdas, Eventos, 9.0, addvecinos=false)
        
    saveXV(Caja, "F")                                                                                        
    if D == 2
        if phi < 0.65 
            tdif = 40.0
        else 
            tdif = 15.0
        end
    elseif D == 3
        tdif = 50.0
    end
                                                                                                
    println("Midiendo Difusion")
    medirnDifusion(Caja, Celdas, Eventos, crecimiento, phi, tdif)
    saveXV(Caja, "Dif")
end   
    
function medirnDifusion(Caja::caja{D}, Celdas::celdas{D}, Eventos::eventos, crecimiento::Float64, phi::Float64, tdif::Real) where D
    DatosNombreArchivo = datosNombreArchivo(Caja.N, Caja.D, Caja.Ntipos, phi, crecimiento)
    dt = 1.0
    NombreArchivoD = "Difdt$(dt)"*DatosNombreArchivo*".csv"

    esferasIniciales = deepcopy(Caja.esferas)
       
    tiempoInicial = Caja.tiempoG
    tiempoFinal = tiempoInicial + tdif
    
    MSDS = Float64[]
    tiempos = Float64[]
    NVMP = Int64[]
    
    prints = 1
    while Caja.tiempoG < tiempoFinal
        if prints*tdif/5 < (Caja.tiempoG - tiempoInicial)
            println("Avance de difusion a $prints / 5")
            prints += 1
        end
        Caja, Celdas, Eventos = evolucionarTconDT(Caja, Celdas, Eventos, 0.1, dt)
        Caja = mover_esferas(Caja)
        
        dt = Caja.tiempoG - tiempoInicial
        MSD = meanSquareDisplacement(esferasIniciales, Caja)
            
        nesferas, nvecinos = histogramaMono(1, Caja)
        m = maximum(nesferas)
        nv = 0
        ancho = 0
        for i in 1:Caja.N
            if nesferas[i] == m
                nv = i
            end
        end
        nvmp = nvecinos[nv]
            
        push!(MSDS, MSD)
        push!(tiempos, dt)
        push!(NVMP, nvmp)
                                                   
    end
    CSV.write(NombreArchivoD, DataFrame(MSD = MSDS, NVMP = NVMP, tiempo = tiempos))
end
    
function evolucionarTconDT(Caja::caja{D}, Celdas::celdas{D}, Eventos::eventos, tiempo::Float64, dt::Float64) where D
    t = Caja.tiempoG 
    while Caja.tiempoG < t + tiempo
        Caja, Celdas, Eventos = proximo_evento(Caja, Celdas, Eventos, agregarvecinos=true, dtV = dt)
    end
    Caja, Celdas, Eventos
end
    
    
function proximo_eventoC(Caja::caja{D}, Celdas::celdas{D}, Eventos::eventos,  contCol::Int64 = 0; vecinos::Bool = true, dtV = 0.0) where D
    tiempoMin = minimum(Eventos.tiemposCyT)
    n1, j = findfirst(Eventos.tiemposCyT, tiempoMin)
    
    Caja.tiempoG = tiempoMin
    n2 = Eventos.eventoCyT[n1, j]

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
            contCol += 1
            
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
    Caja, Celdas, Eventos, contCol
end
            
            
function medirSoloDifusionCol(N::Int64, tipo::Int64, D::Int64, crecimiento::Float64, phi::Float64; propGrCh=1.4)
       # tdif en args
    println("Creciendo particulas con phi = $phi , con $tipo radios en $D D")
    Caja, Celdas, Eventos = CrecerDiscosAPhi(N, tipo, D,  crecimiento, phi, propGC=propGrCh)
    Eventos = eventos(Caja.tiempoG, Caja, Celdas)
    
#     SqrtEf = medirEnergia(Caja, Celdas, Eventos, phi, crecimiento, modNombre="I")
    println("Termalizando")
    Caja, Celdas, Eventos = evolucionarFlu(Caja, Celdas, Eventos, 9.0, addvecinos=false)
        
    saveXV(Caja, "F")      
    
#     SqrtEf = medirEnergia(Caja, Celdas, Eventos, phi, crecimiento, modNombre="F")         
    if D == 2
        
        colDif = 100_000_000
#       colDif = 10*(N^2)
    elseif D == 3
        colDif = 100_000_000
    end
                                                                                                
    println("Midiendo Difusion")
    medirDifusionCol(Caja, Celdas, Eventos, crecimiento, phi, colDif)
    saveXV(Caja, "Dif")
#     SqrtEf = medirEnergia(Caja, Celdas, Eventos, phi, crecimiento, modNombre="Dif")
end   
                
function evolucionarNumeroCol(Caja::caja{D}, Celdas::celdas{D}, Eventos::eventos, nCol::Int64; dt::Float64 = 0.0) where D
    nColAcum = 0
    while nColAcum < nCol
        if dt == 0.0
            Caja, Celdas, Eventos, nColAcum = proximo_eventoC(Caja, Celdas, Eventos, nColAcum)# agregarvecinos=true, dtV = dt
        else
            Caja, Celdas, Eventos, nColAcum = proximo_eventoC(Caja, Celdas, Eventos, nColAcum, agregarvecinos=true, dtV = dt)
        end
    end
    Caja, Celdas, Eventos
end
                
function medirDifusionCol(Caja::caja{D}, Celdas::celdas{D}, Eventos::eventos, crecimiento::Float64, phi::Float64, colDif::Int) where D
    DatosNombreArchivo = datosNombreArchivo(Caja.N, Caja.D, Caja.Ntipos, phi, crecimiento)
    NombreArchivoD = "Dif-nCol100M"*DatosNombreArchivo*".csv"
    NombreArchivoN12 = "N1N2-"*DatosNombreArchivo*".csv"

    esferasIniciales = deepcopy(Caja.esferas)
    colDifAcum = 0
    tiempoInicial = Caja.tiempoG
                        
    MSDS = Float64[]
    numCols = Int64[]
    tiempos = Float64[]
    NVMP = Int64[]
    AvgVs = Float64[]
                            
    AvgVPPs = Float64[]
    AvgVPGs = Float64[]
    AvgVGGs = Float64[]
    
    pasosMSD = 50_000
    first = true
    if Caja.Ntipos == 1         
        N1, N2 = (0, 0)         
    elseif Caja.Ntipos == 2
        N1, N2 = contarGrCh(Caja.esferas)
    end
#     dtV = 0.0
                    
    prints = 1
    while colDifAcum < colDif
        if prints*colDif/5 < colDifAcum
            println("Avance de difusion a $prints / 5")
            prints += 1
        end
        if first
            TC = false
            while TC == false
                Caja, Celdas, Eventos = proximo_evento(Caja, Celdas, Eventos)
                for n1 in 1:Caja.N
                    if Caja.numeroColisiones[n1] < 1
                        TC = false
                        break
                    else
                        TC = true
                    end
                end
            end
            println("Todas colisionaron")
    
            first = false                                
        else      
            Caja, Celdas, Eventos =  evolucionarNumeroCol(Caja, Celdas, Eventos, pasosMSD) #reseteo vecinos dt = dtV    
            colDifAcum += pasosMSD     
        end
                                    
                                    
        Caja = mover_esferas(Caja)
        MSD = meanSquareDisplacement(esferasIniciales, Caja)
            
        AvgVec = sum(Caja.numeroVecinos)/Caja.N
        nesferas, nvecinos = histogramaMono(1, Caja)
        m = maximum(nesferas)
        nv = 0
        ancho = 0
        for i in 1:Caja.N
            if nesferas[i] == m
                nv = i
            end
        end
        nvmp = nvecinos[nv]
        
        AvgVecPP, AvgVecPG, AvgVecGG = (0.0, 0.0, 0.0)
        if Caja.Ntipos == 1
            nothing
        elseif Caja.Ntipos == 2
            for n1 in 1:Caja.N
                if Caja.esferas[n1].tipo == 1
                    AvgVecPP += Caja.tipoDeVecinos[n1, 1]
                    AvgVecPG += Caja.tipoDeVecinos[n1, 2]
                elseif Caja.esferas[n1].tipo == 2
                    AvgVecPG += Caja.tipoDeVecinos[n1, 1]
                    AvgVecGG += Caja.tipoDeVecinos[n1, 2]
                end
            end
            AvgVecPP = sum(AvgVecPP)/N1
            AvgVecPG = sum(AvgVecPG)/(N1+N2)
            AvgVecGG = sum(AvgVecGG)/N2
        end

        push!(MSDS, MSD)
        push!(numCols, colDifAcum)
        push!(tiempos, Caja.tiempoG - tiempoInicial)
        push!(NVMP, nvmp)
        push!(AvgVs, AvgVec)
        push!(AvgVPPs, AvgVecPP)
        push!(AvgVPGs, AvgVecPG)       
        push!(AvgVGGs, AvgVecGG)                                    
    end
    CSV.write(NombreArchivoD, DataFrame(MSD = MSDS, NVMP = NVMP, nCols = numCols, tiempos = tiempos, AvgVecinos = AvgVs, AvgVecinosPP = AvgVPPs, AvgVecinosPG = AvgVPGs, AvgVecinosGG = AvgVGGs))
    CSV.write(NombreArchivoN12, DataFrame(N1 = N1, N2 = N2))                                 
end