arraydearrayint(m::Int64) = [Int[] for i in 1:m, j in 1:m]
arraydearraytup(m::Int64) = [Array{Array{Int64,1},1}() for i in 1:m, j in 1:m]

function celdas_vecinas_iniciales(m::Int64)
    cvi = arraydearraytup(m)
    for i in 1:m
        for j in 1:m
            for h in -1:0, l in -1:1
                if h == 0 && l == -1  else push!(cvi[i,j], [i+h,j+l]) end
            end
        end
    end
    
    for n in 1:m^2
        for v in 1:5
            for i in 1:2
                if cvi[n][v][i] == 0
                    cvi[n][v][i] = m
                end
                if cvi[n][v][i] == m+1
                    cvi[n][v][i] = 1
                end
            end
        end
    end
    cvi
end



function celdas_vecinas(m::Int64)
    cv = arraydearraytup(m)
    for i in 1:m
        for j in 1:m
            for h in -1:1, l in -1:1
                push!(cv[i,j], [i+h,j+l])
            end
        end
    end
    
    for n in 1:m^2
        for v in 1:9
            for i in 1:2
                if cv[n][v][i] == 0
                    cv[n][v][i] = m
                end
                if cv[n][v][i] == m+1
                    cv[n][v][i] = 1
                end
            end
        end
    end
    cv
end

function celdas_frontera(m::Int64)
    cf = Array{Array{Int64,1},1}(undef, 0)
    for i in 1:m
        for j in 1:m
            if i == 1
                push!(cf, [i,j])
            elseif j == 1 && i != 1
                push!(cf, [i,j])
            elseif i == m 
                push!(cf, [i,j])
            elseif j == m  && i != m 
                push!(cf, [i,j])
            end
        end
    end
    cf
end

arraydearrayint3d(m::Int64) = [Int64[] for i in 1:m, j in 1:m, k in 1:m]
arraydearraytup3d(m::Int64) = [Array{Array{Int64,1},1}() for i in 1:m, j in 1:m, k in 1:m]


function celdas_vecinas_iniciales3d(m::Int64)
    cvi = arraydearraytup3d(m)
    for i in 1:m
        for j in 1:m
            for k in 1:m
                for h in -1:0, l in -1:1, g in -1:1
                    if h == 0 && l == -1  
                    else 
                        cvp = [i+h,j+l,k+g]
                        for q in 1:3
                            if cvp[q] == 0
                                cvp[q] = m
                            elseif cvp[q] == m+1
                                cvp[q] = 1
                            end
                        end 
                        push!(cvi[i,j,k], cvp) 
                    end
                end
            end
        end
    end
    cvi
end

function celdas_vecinas3d(m::Int64)
    cv = arraydearraytup3d(m)
    for i in 1:m
        for j in 1:m
            for k in 1:m
                for h in -1:1, l in -1:1, g in -1:1
                    cvp = [i+h,j+l,k+g]
                    for q in 1:3
                        if cvp[q] == 0
                            cvp[q] = m
                        elseif cvp[q] == m+1
                            cvp[q] = 1
                        end 
                    end 
                    push!(cv[i,j,k], cvp)
                end
            end
        end
    end
    cv
end

function celdas_frontera3d(m::Int64)
    cf = Array{Array{Int64,1},1}(undef, 0)
    for i in 1:m
        for j in 1:m
            for k in 1:m
                if i == 1
                    push!(cf, [i,j,k])
                elseif j == 1 && i != 1
                    push!(cf, [i,j,k])
                elseif k == 1 && i != 1
                    push!(cf, [i,j,k])
                elseif i == m 
                    push!(cf, [i,j,k])
                elseif j == m  && i != m
                    push!(cf, [i,j,k])
                elseif k == m && i != m
                    push!(cf, [i,j,k])
                end
            end
        end
    end
    cf
end

import Base: getindex

function getindex(A::Array{Array{Array{Int64,1},1},2}, B::Array{Int64,1})
    A[B[1],B[2]]
end

function getindex(A::Array{Array{Array{Int64,1},1},3}, B::Array{Int64,1})
    A[B[1],B[2],B[3]]
end


function getindex(A::Array{Array{Int64,1},2}, B::Array{Int64,1})
    A[B[1],B[2]]
end

function getindex(A::Array{Array{Int64,1},3}, B::Array{Int64,1})
    A[B[1],B[2],B[3]]
end

function corregirExtremos(A::Array{Int64,1}, m::Int64, x::SVector{2, Float64})
    D = length(A)                         
    correcionExtremo = zeros(D)
    for i in 1:D
        if A[i] == 0
            correcionExtremo[i] += 1.0
            A[i] = m
        elseif A[i] == m+1
            correcionExtremo[i] -= 1.0
            A[i] = 1
        end
    end
    A, x + correcionExtremo
end
                                                
mutable struct celdas{D}
    m::Int64 #numero de celdas por lado sin contar celdas de esferas imagen
    l::Float64
    residentes::Array{Array{Int64,1},D}
    presenteEnCelda::Array{Array{Int64,1},1}
    vecinas::Array{Array{Array{Int64,1},1},D}
    vecinasIniciales::Array{Array{Array{Int64,1},1},D}
    frontera::Array{Array{Int64,1},1}
    LC::Array{Float64,1} #limites celdas x y
    function celdas(Caja::caja{D}) where D
        if D == 2
            if Caja.Ntipos == 1
                l = 2*sqrt(0.906/(Caja.N*pi))
            elseif Caja.Ntipos == 2
                N1, N2 = contarGrCh(Caja.esferas)
                l = 2*sqrt(0.80 / (pi*(N1 + (N2*Caja.propGCrn[1]^2) ) ))
#                 l *= Caja.propGCrn[1]*1.03
                l *= Caja.propGCrn[1]*1.1
            end
        elseif D == 3
            l = 2* (((3/4)*0.74)/(Caja.N*pi))^(1/3)
        end
        M::Int64 = floor(Caja.L/l) # numero de celdas por dimension
        LC = LinRange(0, 1, M+1)
        l = LC[2]-LC[1]
        if D == 2
            Residentes = arraydearrayint(M)
        elseif D == 3
            Residentes = arraydearrayint3d(M)
        end
        PresenteEnCelda = Array{Array{Int64,1},1}(undef, Caja.N)
        
        for n1 in 1:Caja.N
            ijk::Array{Int64,1} = cld.(Caja.esferas[n1].x, l)
            ijk, Caja.esferas[n1].x = corregirExtremos(ijk, M, Caja.esferas[n1].x)
            push!(Residentes[ijk], n1)
            PresenteEnCelda[n1] = ijk
        end
        if D == 2        
            cv = celdas_vecinas(M)
            cvi = celdas_vecinas_iniciales(M)
            cf = celdas_frontera(M)
        elseif D == 3
            cv = celdas_vecinas3d(M)
            cvi = celdas_vecinas_iniciales3d(M)
            cf = celdas_frontera3d(M)
        end
        new{D}(M, l, Residentes, PresenteEnCelda, cv, cvi, cf, LC)
    end
end