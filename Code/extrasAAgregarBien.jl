########### Arreglar
function makeTuple(x::SVector{2, Float64})
    (x[1],x[2])
end

function particulasImagen(Caja::caja{2}, Celdas::celdas{2})
    imagenes = SVector{2, Float64}[]
    for C in Celdas.frontera
        for n1 in Celdas.residentes[C]
            X, x = particulaImagen(n1, Caja, Celdas)
            push!(imagenes, x[1])
            push!(imagenes, x[2])
            push!(imagenes, X)            
        end
    end
    unique(imagenes)
end

function getSitesVoronoi(Caja::caja{2}, Celdas::celdas{2})
    sites = Tuple{Float64, Float64}[]
    NIS = particulasImagen(Caja, Celdas)
    
    for i in 1:Caja.N
        push!(sites, makeTuple(Caja.esferas[i].x))
    end
    for i in 1:length(NIS)
        push!(sites, makeTuple(NIS[i]))
    end
    unique(sites)
end

function getVoronoiNeighbors(sites::Vector{Tuple{Float64, Float64}})
    include("/home/aramis/Voronoi/voronoi.jl")
    voronoi = getVoronoiDiagram(sites)
    
    neighborsnumber = Array{Int64, 1}[]
    for n1 in 1:length(voronoi.faces)
        face = voronoi.faces[n1]
        tmpn = []
        halfedge = face.outerComponent
        while true
            neighbor = halfedge.twin.incidentFace.site

            n2 = 0
            for i in 1:length(sites)
                if sites[i] == neighbor
                    n2 = i
                    break
                end
            end
            if n2 != 0.0
                push!(tmpn, n2)
            end
            
            halfedge = halfedge.next
            if halfedge === face.outerComponent
                break
            end
        end
        push!(neighborsnumber, tmpn)
    end
    neighborsnumber
end

function angulo(n1::Int64, n2::Int64, sites::Vector{Tuple{Float64, Float64}})
    X = sites[n1] .- sites[n2]
#     θ = atan(X[2], X[1])    
    θ = mod(atan(X[2], X[1]), 2pi)
end

function psi6Voronoi(Caja::caja{2}, Celdas::celdas{2}) 
    PO6 = 0
    sites = getSitesVoronoi(Caja, Celdas)
    neighbors = getVoronoiNeighbors(sites)
    
    for n1 in 1:Caja.N
        psi = 0
        for n2 in neighbors[n1]
            psi += exp(im*6*angulo(n1, n2, sites))
        end
        psi /= length(neighbors[n1])
        PO6 += norm(psi)
    end
    PO6 /= Caja.N
    PO6 
end

###### Esto funciona bien

function transforma_vector_tupla(arreglo)
    [(arreglo[i][1], arreglo[i][1]) for i in 1:length(arreglo)]
end

function makeTuple(x::SVector{2, Float64})
    (x[1],x[2])
end

function particulasImagen(Caja::caja{2}, Celdas::celdas{2})
    imagenes = SVector{2, Float64}[]
    for C in Celdas.frontera
        for n1 in Celdas.residentes[C]
            X, x = particulaImagen(n1, Caja, Celdas)
            push!(imagenes, x[1])
            push!(imagenes, x[2])
            push!(imagenes, X)            
        end
    end
    unique(imagenes)
end


function getSitesVoronoi(Caja::caja{2}, Celdas::celdas{2})
    sites = Tuple{Float64, Float64}[]
    NIS = particulasImagen(Caja, Celdas)
    
    for i in 1:Caja.N
        push!(sites, makeTuple(Caja.esferas[i].x))
    end
    for i in 1:length(NIS)
        push!(sites, makeTuple(NIS[i]))
    end
    unique(sites)
end

function elimina_frontera(voronoi)
    n = length(voronoi.faces)
    I = []
    for i in 1:n
        if voronoi.faces[i].area < Inf
            push!(I, i)
        end
    end
    return I
end

function elimina_centro(voronoi)
    n = length(voronoi.faces)
    I = []
    for i in 1:n
        if voronoi.faces[i].area == Inf
            push!(I, i)
        end
    end
    return I
end

function vecinos(voronoi, i)
    vecinosi = []
    face = voronoi.faces[i]
    halfedge = face.outerComponent
    while true
        p0 = halfedge.origin.coordinates
        p1 = halfedge.twin.origin.coordinates
        neighbor = halfedge.twin.incidentFace.site
        push!(vecinosi, neighbor)
        halfedge = halfedge.next
        if halfedge === face.outerComponent
            break
        end
    end
    return vecinosi
end
function Ψ_local(n, xs, vecinosi)
    xx = [[vecinosi[j][1]-xs[1],vecinosi[j][2]-xs[2]]  for j in 1:length(vecinosi)]
    θ_i_j = [mod(atan(xx[j][2], xx[j][1]), 2π) for j in 1:length(xx)]
    norm(sum(exp(im*n*θ_i_j[i]) for i in 1:length(xx)))/length(xx)
end

function Ψ_local2(n, xs, vecinosi)
    xx = [[vecinosi[j][1]-xs[1],vecinosi[j][2]-xs[2]]  for j in 1:length(vecinosi)]
    θ_i_j = [mod(atan(xx[j][2], xx[j][1]), 2π) for j in 1:length(xx)]
    sum(exp(im*n*θ_i_j[i]) for i in 1:length(xx))/length(xx)
end

function obtén_vecinos(sitios)
    voronoi = getVoronoiDiagram(sitios);
    I = elimina_frontera(voronoi)
    vecinos_i = []
    for i in 1:length(sitios)
        vs = vecinos(voronoi, i)
        push!(vecinos_i, vs)
    end
    return I, vecinos_i, [voronoi.faces[i].site for i in 1:length(sitios)]
end
   
    
function Ψ(n, sitios)
    N = length(sitios)
    I, vecinos, sites = obtén_vecinos(sitios)
    sum(Ψ_local(n, sites[i], vecinos[i]) for i in I)/length(I) 
end

function Ψ_global(n, sitios)
    N = length(sitios)
    I, vecinos, sites = obtén_vecinos(sitios)
    norm(sum(Ψ_local2(n, sites[i], vecinos[i]) for i in I))/length(I) 
end

function caja(xrt::DataFrame, v::DataFrame, datos::DataFrame)
    L = datos.cajon[2] - datos.cajon[1]
    N = datos.N[1]
    D = datos.D[1]
    Esferas = Array{esfera{D}, 1}(undef, N)
    vecinos = Array{Array{vecinosDinamicos,1},1}(undef, N)
    transferenciasPeriodicas = Array{Array{Int64,1},1}(undef, N)
    numeroVecinos = zeros(Int64, N)
    numeroColisiones = zeros(Int64, N)
    tipoDeVecinos = zeros(Int64, N, 2)
    crecimientos = [0.0]
    
    for n1 in 1:N
        x = SVector{D,Float64}([xrt.x1[n1], xrt.x2[n1]])
        vel = SVector{D,Float64}([v.v1[n1], v.v2[n1]])
        if datos.Ntipos[1] == 1
            Sn1 = esfera{D}(x, vel, xrt.r[n1], 1.0, 0.0, 0.0, 0)
        elseif datos.Ntipos[1] == 2
            masa = 1.0
            if xrt.tipo == 2
                try
                    masa = (datos.propGrCh[1])^2
                catch
                    masa = (datos.propGCrn[1])^2
                end
            end
            Sn1 = esfera{D}(x, vel, xrt.r[n1], masa, 0.0, 0.0, xrt.tipo)
        elseif datos.Ntipos[1] > 2
                error("Aun nada")
        end
        Esferas[n1] = Sn1
        vecinos[n1] = []
        transferenciasPeriodicas[n1] = zeros(Int, D)
    end
                
    prop = [1.0, 1.0]
    try
        prop = [datos.propGrCh[1], 2.0]
    catch
        prop = [datos.propGCrn[1], datos.propGCrn[2]]
    end
        
 
    caja{D}(N, datos.Ntipos[1], D, Esferas, 0.0, prop, crecimientos, vecinos, numeroVecinos, numeroColisiones, tipoDeVecinos, transferenciasPeriodicas, [datos.cajon[1], datos.cajon[2]], L)
end
                