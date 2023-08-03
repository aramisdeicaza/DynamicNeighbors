struct evento
    tiempo::Float64
    n1::Int64
    n2::Int64
    T::Int64
end

import Base: +, -, first, isless, findfirst

function isless(a::evento, b::evento)
    return a.tiempo < b.tiempo
end

#probar la de julia otra vez para diferencia de conjuntos
function -(A::Array{Array{Int64,1},1}, B::Array{Array{Int64,1},1})
    diferencia = Array{Int64,1}[]
    for a in A
        if a ∉ B
            push!(diferencia, a)
        end
    end
    diferencia
end

function +(a::evento, b::Float64)
    evento(a.tiempo+b, a.n1, a.n2, a.T)
end

function tiempo_colision(n1::Int64, n2::Int64, Caja::caja{D}, Celdas::celdas{D}) where D
    s1 = Caja.esferas[n1]
    s2 = Caja.esferas[n2]
    
    p1 = s1.x + s1.v.*(Caja.tiempoG - s1.t)
    p2 = s2.x + s2.v.*(Caja.tiempoG - s2.t)
    Δp = p1 - p2
    Δp = distanciaPorComponente(Δp, Caja) 
    
    v1 = s1.v
    v2 = s2.v
 
    Δv = v1 - v2
    
    σ = s1.c + s2.c

    if σ == 0.0
        ρ = s1.r + s2.r
    else
        # ρ = Caja.tiempoG*σ 
        ρ = s1.r + s2.r + s1.c*(Caja.tiempoG - s1.t) + s2.c*(Caja.tiempoG - s2.t)
    end
    a, b, c = (norm(Δv)^2 - σ^2, Δv⋅Δp - σ*ρ, norm(Δp)^2 - ρ^2)
    if c < -1e-10*ρ #funciona con -1e-6
        if 2*s1.r > Celdas.l || 2*s2.r > Celdas.l
            error("Alguna esfera ha crecido mas que la celda")
        end
        # println("Traslape encontrado")
#         @show n1, n2, a, b, c
#         error("Traslape encontrado, $c")
    end
    t = quadratic(a, b, c)
end

function quadratic(a::Float64, b::Float64, c::Float64)
    d2 = b^2 - a*c
    t = Inf
    if c <= 0.0 && b < 0.0
        t = 0.0
    elseif d2 >= -10*eps()
        if (d2 < 0.0); d2 = 0.0; end
        if b < 0.0
            t = c / (-b + sqrt(d2))
        elseif a < 0.0 && b > 0.0
            t = -(b + sqrt(d2))/a
        end
    else 
        t = Inf
    end
    t
end

function tiempo_pared(n1::Int64, Caja::caja{D}) where D
    prov = BinaryMinHeap{evento}()
    
    s1 = Caja.esferas[n1]
    pn1 = s1.x + s1.v.*(Caja.tiempoG - s1.t)
    rnow = s1.r + s1.c*(Caja.tiempoG - s1.t)
    
    for i in 1:D # -1 derecha, -2 izquierda, -3 arriba, -4 abajo, -5 zarriba, -6 zabajo
        v1 = s1.v[i] + s1.c
        v2 = s1.v[i] - s1.c
        if v1 > 0.0
            pared = -1 - 2*(i-1)
            tiempo = (Caja.caja[2] - rnow - pn1[i]  ) / v1
        elseif v2 < 0.0 
            pared = -2 - 2*(i-1)
            tiempo = (Caja.caja[1] + rnow - pn1[i] ) / v2
        else
            pared = 0
            tiempo = Inf
        end
        push!(prov, evento(tiempo, n1, pared, 0))
    end
    return top(prov)
end

function tiempos_iniciales_colision(Caja::caja{D}, Celdas::celdas{D}; pared=false) where D
    ts = zeros(Caja.N)
    vs = zeros(Int64, Caja.N)
    
    for n1 in 1:Caja.N
        prov = BinaryMinHeap{evento}()
        ijk = Celdas.presenteEnCelda[n1]
        
        #checa Caja en las celdas vecinas iniciales
        for C in Celdas.vecinasIniciales[ijk]
            for n2 in Celdas.residentes[C]
                if n1 != n2
                    tiempo = tiempo_colision(n1, n2, Caja, Celdas)
                else
                    tiempo = Inf
                end
                push!(prov, evento(tiempo, n1, n2, 0))
            end
        end
        
        #checa colision con paredes si esta presente en celda de frontera
        if pared && ijk ∈ Celdas.frontera
            push!(prov, tiempo_pared(n1, Caja))
        end
        
        e = top(prov)
        ts[n1] = e.tiempo
        vs[n1] = e.n2
    end
    ts, vs
end

function tiempo_transferencia(n1::Int64, Caja::caja{D}, Celdas::celdas{D}) where D
    ijk = Celdas.presenteEnCelda[n1]
    prov = BinaryMinHeap{evento}()
    
    s1 = Caja.esferas[n1]
    pn1 = s1.x + s1.v.*(Caja.tiempoG - s1.t)
    
    for i in 1:D # 1 derecha, 2 izquierda, 3 arriba, 4 abajo
        if s1.v[i] > 0.0
            transferencia = 1 + 2*(i-1)
            tiempo = (Celdas.LC[ijk[i]+1] - pn1[i]) / s1.v[i]
        elseif s1.v[i] < 0.0 
            transferencia = 2 + 2*(i-1)
            tiempo = (Celdas.LC[ijk[i]] - pn1[i]) / s1.v[i]
        else
            transferencia = 0
            tiempo = Inf
        end
        push!(prov, evento(tiempo, n1, 0, transferencia))
    end
    return top(prov)
end

function tiempos_iniciales_transferencia(Caja::caja{D}, Celdas::celdas{D}) where D
    ts = zeros(Caja.N)
    transferencia = zeros(Int64, Caja.N)
    
    for n1 in 1:Caja.N
        e = tiempo_transferencia(n1, Caja, Celdas)
        ts[n1] = e.tiempo
        transferencia[n1] = e.T
    end
    ts, transferencia
end

mutable struct eventos
    tiemposCyT::Array{Float64, 2}
    eventoCyT::Array{Int64, 2}
    # function eventos(nombreDir::String)
    #     df = leerEventos(nombreDir)
    #     N = length(df.tiemposC)
    #     ts = zeros(N, 2)
    #     ev = zeros(Int64, N, 2)
    #     ts[:,1], ev[:,1] = df.tiemposC, df.eventosC
    #     ts[:,2], ev[:,2] = df.tiemposT, df.eventosT
    #     new(ts, ev)
    # end
end
                                    
function eventos(Caja::caja{D}, Celdas::celdas{D}) where D
    ts = zeros(Caja.N, 2)
    ev = zeros(Int64, Caja.N, 2)
    ts[:,1], ev[:,1] = tiempos_iniciales_colision(Caja, Celdas)
    ts[:,2], ev[:,2] = tiempos_iniciales_transferencia(Caja, Celdas)
    eventos(ts, ev)
end         
                                                        
function eventos(t0::Float64, Caja::caja{D}, Celdas::celdas{D}) where D
    # Para usar depuedes de normalizar velocidades o parar o dar crecimiento
    ts = zeros(Caja.N, 2)
    ev = zeros(Int64, Caja.N, 2)
    ts[:,1], ev[:,1] = tiempos_iniciales_colision(Caja, Celdas)
    ts[:,2], ev[:,2] = tiempos_iniciales_transferencia(Caja, Celdas)
    ts = ts .+ t0
    eventos(ts, ev)
end

function tiempo_colision_PTransferencia(n1::Int64, Celdaold::Array{Int64,1}, Caja::caja{D}, Celdas::celdas{D}, Eventos::eventos; pared=false) where D
    prov = BinaryMinHeap{evento}()
    push!(prov, evento(Inf, n1, n1, 0))
    
    ijk = Celdas.presenteEnCelda[n1]
    nuevasVecinas = Celdas.vecinas[ijk] - Celdas.vecinas[Celdaold]
    
    for C in nuevasVecinas
        for n2 in Celdas.residentes[C]
            if n1 != n2
                tiempo = tiempo_colision(n1, n2, Caja, Celdas)
            else 
                tiempo = Inf
            end
            push!(prov, evento(tiempo, n1, n2, 0))
        end
    end

    if pared && ijk ∈ Celdas.frontera
        push!(prov, tiempo_pared(n1, Caja))
    end  
    
    e = top(prov) + Caja.tiempoG
    
    if Eventos.tiemposCyT[n1, 1] > e.tiempo
        Eventos.tiemposCyT[n1, 1] = e.tiempo
        Eventos.eventoCyT[n1, 1] = e.n2
    end

    Eventos
end

function tiempos_colision_nuevos(n1s, Caja::caja{D}, Celdas::celdas{D}, Eventos::eventos; pared = false) where D
    for n1 in n1s
        prov = BinaryMinHeap{evento}()
        ijk = Celdas.presenteEnCelda[n1]
        
        for C in Celdas.vecinas[ijk]
            for n2 in Celdas.residentes[C]
                if n2 ∉ n1s
                    tiempo = tiempo_colision(n1, n2, Caja, Celdas)
                else 
                    tiempo = Inf
                end
                push!(prov, evento(tiempo, n1, n2, 0))
            end
        end

        if pared && ijk ∈ Celdas.frontera
            push!(prov, tiempo_pared(n1, Caja))
        end
        
        e = top(prov) + Caja.tiempoG

        Eventos.tiemposCyT[n1, 1] = e.tiempo
        Eventos.eventoCyT[n1, 1] = e.n2
    end
    
    Eventos
end

function tiempo_transferencia_nuevo(n1::Int64, Caja::caja{D}, Celdas::celdas{D}, Eventos::eventos) where D
    e = tiempo_transferencia(n1, Caja, Celdas) + Caja.tiempoG
    
    Eventos.tiemposCyT[n1, 2] = e.tiempo
    Eventos.eventoCyT[n1, 2] = e.T
    
    Eventos
end

function recalcular_tiempo_colision(n1s, Caja::caja{D}, Celdas::celdas{D}, Eventos::eventos) where D
    for n1 in n1s
        ijk = Celdas.presenteEnCelda[n1]
        for C in Celdas.vecinas[ijk]
            for n2 in Celdas.residentes[C]
                if Eventos.eventoCyT[n2,1] ∈ n1s
                    Eventos = tiempos_colision_nuevos(n2, Caja, Celdas, Eventos)
                end
            end
        end
    end 
    Eventos
end
