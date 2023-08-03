function circulos(Caja::caja{2})
    θ = 0:0.05*pi:2pi
    nθ = length(θ)
    circulo = zeros(nθ, 2)

    if Caja.Ntipos == 1
        P = 0
    else
        P = 1 / (1 + Caja.propGCrn[2])
    end
    #[t1 = RGB(0.4, 0.55, 0.3), t2 RGB(0.4, 0.55, 0.3+0.33)], 
   # generador de círculos alrededor de la partícula
   for n1 in 1:Caja.N
       circulo[1:nθ, 1] = collect(Caja.esferas[n1].x[1] .+ Caja.esferas[n1].r*cos.(θ'))'
       circulo[1:nθ, 2] = collect(Caja.esferas[n1].x[2] .+ Caja.esferas[n1].r*sin.(θ'))'
        if Caja.esferas[n1].tipo == 2
            plot!(circulo[:,1],circulo[:,2], c=:black, fillrange=0, fillcolor = RGB(0.4, 0.45, 0.3+P))
        else
            plot!(circulo[:,1],circulo[:,2], c=:black, fillrange=0, fillcolor = RGB(0.4, 0.45, 0.3))
        end
    end
    plot!()
end


function circulos3D(Caja::caja{3})
    Num = 32
    u = LinRange(0, 2π, Num)
    v = LinRange(0, π, Num)
    
    for n1 in 1:Caja.N
        s = Caja.esferas[n1]
        x = s.x[1] .+ s.r*cos.(u)*sin.(v)'
        y = s.x[2] .+ s.r*sin.(u)*sin.(v)'
        z = s.x[3] .+ s.r*ones(Num)*cos.(v)'

        plot!(x,y,z, st = :surface,  fz=z, colorbar=false)
    end
    plot!()
end

function particulaImagen(n1::Int64, Caja::caja{2}, Celdas::celdas{2})
    X = Caja.esferas[n1].x
    correccion = zeros(2)
    for i in 1:2
        if X[i] > Celdas.LC[end-1]
            correccion[i] = - Caja.L
        elseif X[i] < Celdas.LC[2]
            correccion[i] =  Caja.L           
        end
    end
    X1 = X + correccion
    
    correcionesquinas = [zeros(2) for i in 1:2]
    if Celdas.presenteEnCelda[n1] ∈ [[1,1], [1,Celdas.m], [Celdas.m,1], [Celdas.m,Celdas.m]]
        for i in 1:2
            if X[i] > Celdas.LC[end-1]
                correcionesquinas[i][i] = - Caja.L
            elseif X[i] < Celdas.LC[2]
                correcionesquinas[i][i] = Caja.L               
            end
        end
    end
    X1, [X + correcionesquinas[1], X + correcionesquinas[2]]
end

function circulosImagen(Caja::caja{2}, Celdas::celdas{2})
    #rdp = radio particulas
    θ = 0:0.05*pi:2pi
    nθ = length(θ)
    circulo = zeros(nθ, 2)
    
    if Caja.Ntipos == 1
        P = 0
    else
        P = 1 / (1 + Caja.propGCrn[2])
    end
#     [t1 RGB(0.4, 0.45, 0.4), t2 = RGB(0.4+P/6, 0.45+P/6, 0.4+2P/3)]
    
    for C in Celdas.frontera
        for n1 in Celdas.residentes[C]
            X, x = particulaImagen(n1, Caja, Celdas)
            if C ∈ [[1,1], [1,Celdas.m], [Celdas.m,1], [Celdas.m,Celdas.m]]
                for i in 1:2
                    circulo[1:nθ, 1] = collect(x[i][1] .+ Caja.esferas[n1].r*cos.(θ'))'
                    circulo[1:nθ, 2] = collect(x[i][2] .+ Caja.esferas[n1].r*sin.(θ'))'
                    if Caja.esferas[n1].tipo == 2
                        plot!(circulo[:,1],circulo[:,2], c=:black, fillrange=0, fillcolor = RGB(0.4+P/6, 0.45+P/6, 0.4+2P/3))
                    else
                        plot!(circulo[:,1],circulo[:,2], c=:black, fillrange=0, fillcolor = RGB(0.4, 0.45, 0.4))
                    end
                end
            end
            circulo[1:nθ, 1] = collect(X[1] .+ Caja.esferas[n1].r*cos.(θ'))'
            circulo[1:nθ, 2] = collect(X[2] .+ Caja.esferas[n1].r*sin.(θ'))'
            if Caja.esferas[n1].tipo == 2
                plot!(circulo[:,1],circulo[:,2], c=:black, fillrange=0, fillcolor = RGB(0.4+P/6, 0.45+P/6, 0.4+2P/3))
            else
                plot!(circulo[:,1],circulo[:,2], c=:black, fillrange=0, fillcolor = RGB(0.4, 0.45, 0.4))
            end
        end
    end
    plot!()
end

function grafica(Caja::caja{2}, Celdas::celdas{2})
    Lc = Celdas.l
    l = Caja.L
    phi = round(packingFraction(Caja), digits = 3)
    t = round(Caja.tiempoG, digits = 4)
    gr()
    psi = round(psi6(Caja), digits = 3)

    plot(title = "Time = $t s, phi = $phi , psi6 = $psi", 
        xlims = (-Lc, Caja.L+Lc), ylims = (-Lc, Caja.L+Lc), xformatter=_->"", yformatter=_->"",
    # xticks = -Lc:Lc:(Caja.L+Lc), yticks = -Lc:Lc:(Caja.L+Lc), 
        aspect_ratio=:equal, leg=false)

    circulos(Caja)
    circulosImagen(Caja, Celdas)

    plot!([(0,Caja.L),(Caja.L,Caja.L),(Caja.L,0),(0,0),(0,Caja.L)], c=:orange)
end

function grafica(Caja::caja{3}, Celdas::celdas{3})
    Lc = Celdas.l
    l = Caja.L
    phi = round(packingFraction(Caja), digits = 3)
    t = round(Caja.tiempoG, digits = 4)
#     pyplot()
    q6 = round(Q6(Caja), digits = 3)

    plot(title = "Time = $t s, phi = $phi , Q6 = $q6", 
        xlims = (-0.02, l + 0.02), ylims = (-0.02, l + 0.02), zlim = (-0.02, l + 0.02),# show =:ijulia,
        # xticks=0.0:Lc:1.0, yticks=0.0:Lc:1.0, zticks=0.0:Lc:1.0,
    ylabel = "y", zlabel = "z", xlabel="x", aspect_ratio=:equal, leg=false)

    circulos3D(Esferas)

    plot!([0,0,0],[l,0,0],[l,l,0], c=:red)
    plot!([l,l,0],[l,l,l],[0,l,l], c=:red)
    plot!([0,l,l],[0,0,l],[0,0,0], c=:red)

#     esquinas = [0  0  0;  0  0  1;  0  1  1;  1  1 1;  1  0  1;  1  0  0; 0 0 0; 0 1 0; 1 1 0; 1 1 1; 0 1 1; 0 1 0; 1 1 0; 1 0 0; 1 0 1; 0 0 1]
#     plot!(esquinas[:,1], esquinas[:,2], esquinas[:,3], c=:red)
    
end

function graficaST(Caja::caja{2}, Celdas::celdas{2})
    Lc = Celdas.l
    l = Caja.L
    phi = round(packingFraction(Caja), digits = 4)
    gr()
    psi = round(psi6(Caja), digits = 3)

    plot(title = "phi = $phi , psi6 = $psi", 
        xlims = (-Lc, Caja.L+Lc), ylims = (-Lc, Caja.L+Lc), 
#     xticks = -Lc:Lc:(Caja.L+Lc), yticks = -Lc:Lc:(Caja.L+Lc), 
        aspect_ratio=:equal, leg=false)

    circulos(Caja)
    circulosImagen(Caja, Celdas)
    plot!([(0,Caja.L),(Caja.L,Caja.L),(Caja.L,0),(0,0),(0,Caja.L)], c=:orange)
end
    
function graficaST(Caja::caja{3}, Celdas::celdas{3})
    Lc = Celdas.l
    l = Caja.L
    phi = round(packingFraction(Caja), digits = 3)
#     pyplot()
    q6 = round(Q6(Caja), digits = 3)

    plot(title = "phi = $phi , Q6 = $q6", 
        xlims = (-0.02, l + 0.02), ylims = (-0.02, l + 0.02), zlim = (-0.02, l + 0.02),# show =:ijulia,
#         xticks=0.0:Lc:1.0, yticks=0.0:Lc:1.0, zticks=0.0:Lc:1.0,
    ylabel = "y", zlabel = "z", xlabel="x", leg=false)

    circulos3D(Caja)

    plot!([0,0,0],[l,0,0],[l,l,0], c=:red)
    plot!([l,l,0],[l,l,l],[0,l,l], c=:red)
    plot!([0,l,l],[0,0,l],[0,0,0], c=:red)
        
#     esquinas = [0  0  0;  0  0  1;  0  1  1;  1  1 1;  1  0  1;  1  0  0; 0 0 0; 0 1 0; 1 1 0; 1 1 1; 0 1 1; 0 1 0; 1 1 0; 1 0 0; 1 0 1; 0 0 1]
#     plot!(esquinas[:,1], esquinas[:,2], esquinas[:,3], c=:red)
end
                            
                            
function graficaSP(Caja::caja{2}, Celdas::celdas{2})
    Lc = Celdas.l
    l = Caja.L
    phi = round(packingFraction(Caja), digits = 4)
    gr()

    plot(title = "phi = $phi", 
        xlims = (-Lc, Caja.L+Lc), ylims = (-Lc, Caja.L+Lc), 
#     xticks = -Lc:Lc:(Caja.L+Lc), yticks = -Lc:Lc:(Caja.L+Lc), 
        aspect_ratio=:equal, leg=false)

    circulos(Caja)
    circulosImagen(Caja, Celdas)
    plot!([(0,Caja.L),(Caja.L,Caja.L),(Caja.L,0),(0,0),(0,Caja.L)], c=:orange)
end

function circulosVecinos(Caja::caja{2})
    θ = 0:0.05*pi:2pi
    nθ = length(θ)
    circulo = zeros(nθ, 2)
    gradiente = false
    if sum(Caja.numeroVecinos)/Caja.N > 14
        clims = extrema([minimum(Caja.numeroVecinos), maximum(Caja.numeroVecinos)])
        gradiente = true
        cMin::Int = minimum(Caja.numeroVecinos)
        nColores::Int = maximum(Caja.numeroVecinos) - cMin + 1
    else
        clims = extrema([5, 14])
    end
    
   for n1 in 1:Caja.N
       circulo[1:nθ, 1] = collect(Caja.esferas[n1].x[1] .+ Caja.esferas[n1].r*cos.(θ'))'
       circulo[1:nθ, 2] = collect(Caja.esferas[n1].x[2] .+ Caja.esferas[n1].r*sin.(θ'))'
        n = Caja.numeroVecinos[n1]
        if n > 14
            n = 14
        elseif n < 5
            n = 5
        end
        if gradiente 
            plot!(circulo[:,1],circulo[:,2], c=:black, fillrange=0, 
                fill=(true,palette(:Spectral_10, nColores)[Caja.numeroVecinos[n1]-cMin+1]), zcolor=Caja.numeroVecinos[n1]-cMin+1, clims=clims)
        else
            plot!(circulo[:,1],circulo[:,2], c=:black, fillrange=0, fill=(true,palette(:Spectral_10)[n-4]), zcolor=n-4, clims=clims)
        end
        
    end
    plot!()
end
                                           
function graficaVecinos(Caja::caja{2}, Celdas::celdas{2})
    Lc = Celdas.l
    l = Caja.L
    phi = round(packingFraction(Caja), digits = 3)
    psi = round(psi6(Caja), digits = 3)
    gr()

    plot(title = "phi = $phi , psi6 = $psi", 
        xlims = (-Lc, Caja.L+Lc), ylims = (-Lc, Caja.L+Lc), 
        aspect_ratio=:equal, leg=false)

    gradiente = false
    if sum(Caja.numeroVecinos)/Caja.N > 14
        clims = extrema([minimum(Caja.numeroVecinos), maximum(Caja.numeroVecinos)])
        gradiente = true
    else
        clims = extrema([5, 14])
    end
    
    circulosVecinos(Caja)
    circulosImagen(Caja, Celdas)
    
    plot!([(0,Caja.L),(Caja.L,Caja.L),(Caja.L,0),(0,0),(0,Caja.L)], c=:orange, leg=false)
    if gradiente
        scatter!([-1], [-1], clims=clims, color = cgrad(:Spectral_10), 
            zcolor=[minimum(Caja.numeroVecinos), maximum(Caja.numeroVecinos)], colorbar=true)
    else
        scatter!([-1], [-1], clims=clims, color = palette(:Spectral_10), zcolor=[5, 14], colorbar=true)
    end
    
end
                                                
                                                                                       
function graficaVecinos(Caja::caja{2}, Celdas::celdas{2}, tiempo::Float64)
    Lc = Celdas.l
    l = Caja.L
    phi = round(packingFraction(Caja), digits = 3)
    psi = round(psi6(Caja), digits = 3)
    gr()

    plot(title = "phi = $(phi), t = $(tiempo), psi6 = $(psi)", 
        xlims = (-Lc, Caja.L+Lc), ylims = (-Lc, Caja.L+Lc), 
        aspect_ratio=:equal, leg=false)

    gradiente = false
    if sum(Caja.numeroVecinos)/Caja.N > 14
        clims = extrema([minimum(Caja.numeroVecinos), maximum(Caja.numeroVecinos)])
        gradiente = true
    else
        clims = extrema([5, 14])
    end
    
    circulosVecinos(Caja)
    circulosImagen(Caja, Celdas)
    
    plot!([(0,Caja.L),(Caja.L,Caja.L),(Caja.L,0),(0,0),(0,Caja.L)], c=:orange, leg=false)
    if gradiente
        scatter!([-1], [-1], clims=clims, color = cgrad(:Spectral_10), 
            zcolor=[minimum(Caja.numeroVecinos), maximum(Caja.numeroVecinos)], colorbar=true)
    else
        scatter!([-1], [-1], clims=clims, color = palette(:Spectral_10), zcolor=[5, 14], colorbar=true)
    end
    
end

function graficaVecinos(Caja::caja{2}, Celdas::celdas{2}, t::Float64, dt::Float64, dtV::Float64)
    Lc = Celdas.l
    l = Caja.L
    phi = round(packingFraction(Caja), digits = 3)
    psi = round(psi6(Caja), digits = 3)
    gr()

    plot(title = "t = $t s, dt = $dt s, phi = $phi , psi6 = $psi , dtV = $dtV s", 
        xlims = (-Lc, Caja.L+Lc), ylims = (-Lc, Caja.L+Lc), 
        aspect_ratio=:equal, leg=false)

    gradiente = false
    if sum(Caja.numeroVecinos)/Caja.N > 14
        clims = extrema([minimum(Caja.numeroVecinos), maximum(Caja.numeroVecinos)])
        gradiente = true
    else
        clims = extrema([5, 14])
    end
    
    circulosVecinos(Caja)
    circulosImagen(Caja, Celdas)
    
    plot!([(0,Caja.L),(Caja.L,Caja.L),(Caja.L,0),(0,0),(0,Caja.L)], c=:orange, leg=false)
    if gradiente
        scatter!([-1], [-1], clims=clims, color = cgrad(:Spectral_10), 
            zcolor=[minimum(Caja.numeroVecinos), maximum(Caja.numeroVecinos)], colorbar=true)
    else
        scatter!([-1], [-1], clims=clims, color = palette(:Spectral_10), zcolor=[5, 14], colorbar=true)
    end
    
end

function circulosmucz(Caja::caja{2})
    θ = 0:0.05*pi:2pi
    nθ = length(θ)
    circulo = zeros(nθ, 2)

    clims = extrema([0, 2_000])
    nColores::Int = 2_000
       
   for n1 in 1:Caja.N
        if Caja.numeroVecinos[n1] != 0
            n::Int64 = ceil(Caja.numeroColisiones[n1]/Caja.numeroVecinos[n1])
        else
            n = 1
        end
        circulo[1:nθ, 1] = collect(Caja.esferas[n1].x[1] .+ Caja.esferas[n1].r*cos.(θ'))'
        circulo[1:nθ, 2] = collect(Caja.esferas[n1].x[2] .+ Caja.esferas[n1].r*sin.(θ'))'

        plot!(circulo[:,1],circulo[:,2], c=:black, fillrange=0, 
                fill=(true,palette(:Spectral_10, nColores)[n]), zcolor=n, clims=clims)
        
    end
    plot!()
end

function graficamucz(Caja::caja{2}, Celdas::celdas{2})
    Lc = Celdas.l
    l = Caja.L
    phi = round(packingFraction(Caja), digits = 3)
    t = round(Caja.tiempoG, digits = 4)
    psi = round(psi6(Caja), digits = 3)
    gr()
    plot(title = "phi = $phi , psi6 = $psi", 
        xlims = (-Lc, Caja.L+Lc), ylims = (-Lc, Caja.L+Lc), 
        aspect_ratio=:equal, leg=false)
        
    # Caja.numeroVecinos[n1]/Caja.numeroColisiones[n1]
    clims = extrema([0, 2_000])
    
    circulosmucz(Caja)
    circulosImagen(Caja, Celdas)
    
    plot!([(0,Caja.L),(Caja.L,Caja.L),(Caja.L,0),(0,0),(0,Caja.L)], c=:orange, leg=false)
    scatter!([-1], [-1], clims=clims, color = cgrad(:Spectral_10), 
            zcolor=[0, 2_000], colorbar=true)
end

function animarCrecimientoAPhi(N::Int64, tipo::Int64, D::Int64, crecimiento::Float64, ϕ::Float64, dir::String, frames::Int64; propGrCh::Float64=1.4, v::Float64=5.0)
    Caja, Celdas, Eventos = condiciones_iniciales(N, tipo, D, crecimiento, prop=propGrCh, vm = v)
    tPhi = 0
    if D == 2
        if tipo == 1
            tPhi = (ϕ / (N*pi*crecimiento^2))^(1/2)
        elseif tipo == 2
            N1, N2 = contarGrCh(Caja.esferas)
            tPhi = (ϕ / (pi*N1*crecimiento^2 + pi*N2*(Caja.propGrCh*crecimiento)^2))^(1/2)
        end
    elseif D == 3
        tPhi = ( (3/4)*ϕ / (N*pi*crecimiento^3))^(1/3)
    end
    
    mkdir(dir)
    dt = tPhi/ (frames)
    frame = 1
    
    while Caja.tiempoG < tPhi
        if Caja.tiempoG >= (frame-1)*dt
            Caja = mover_esferas(Caja)
            grafica(Caja, Celdas)
            savefig(dir*"/Fig"*lpad(frame,6,"0")*".png")
            frame += 1
        end
        Caja, Celdas, Eventos = proximo_evento(Caja, Celdas, Eventos, agregarvecinos=false)
    end
#     ffmpeg -framerate 60 -i Fig_%05d.jpg output.mp4
end

function animarDinamicaEnPhi(N::Int64, tipo::Int64, D::Int64, crecimiento::Float64, ϕ::Float64, dtVe::Float64, dir::String, frames::Int64, vecinos::Bool; propGrCh::Float64=1.4, v::Float64=5.0)
    Caja, Celdas, Eventos = CrecerDiscosAPhi(N, tipo, D,  crecimiento, ϕ, propGC=propGrCh)
                                                                        
#     Caja, vm = normalizar_velocidades(Caja)
#     Eventos = eventos(Caja.tiempoG, Caja, Celdas)
    
    dt = (maximum(filter(!isinf,Eventos.tiemposCyT[1,:])) - minimum(Eventos.tiemposCyT) ) / 8
    @show dt
#     dt = round(dt, digits=5)
    frame = 1
    t0 = Caja.tiempoG
    tFinal = Caja.tiempoG + dt*frames
    
    mkpath(dir)
    @show t0, tFinal, Caja.tiempoG                                                
    
    while Caja.tiempoG < tFinal
        if Caja.tiempoG >= t0 + (frame-1)*dt
            Caja = mover_esferas(Caja)
            if D == 2 && vecinos == true
                graficaVecinos(Caja, Celdas, Caja.tiempoG-t0, dt, dtVe)
            else
                grafica(Caja, Celdas)
            end
            savefig(dir*"/Fig"*lpad(frame,6,"0")*".png")
            frame += 1
        end
        Caja, Celdas, Eventos = proximo_evento(Caja, Celdas, Eventos, agregarvecinos=vecinos, dtV = dtVe)
    end
#     ffmpeg -framerate 20 -i Fig%06d.png output.mp4
end