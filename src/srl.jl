#-----------------------------------------------------------------------
module cereal      # cereal module
#-----------------------------------------------------------------------

using LinearAlgebra
using Combinatorics

const RealVec{T<:Real} = Array{T,1}      # Defining vector datatype
const RealMtx{T<:Real} = Array{T,2}      # Defining matrix datatype

#-----------------------------------------------------------------------
function η( V1::RealVec , V2::RealVec )  # Computes Minkowski product
    nv1 = length(V1)
    nv2 = length(V2)

    met = -V1[1]*V2[1]

    if nv1==nv2
        for i=2:nv1
            met += V1[i]*V2[i]
        end
        return met
    else
        print("Vectors are of a different dimension.")
    end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
function ϵ(i,j,k,l)
    return levicivita([i,j,k,l])
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
function ημν( tpfl::DataType=Float64 , dim::Int=4 )
    gη = one(tpfl)*(I(dim))

    gη[1,1] = -gη[1,1]

    return gη
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
function Hodge2( U::RealVec , V::RealVec )
    tpfl=typeof(U[1]) 

    W = zeros(tpfl,4,4)
    for i=1:4
    for j=1:4
        Z  = zero(tpfl)
        for k=1:4
        for l=1:4
            Z += -ϵ(i,j,k,l)*U[k]*V[l]
        end
        end
        W[i,j] = Z
    end
    end

    return W
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
function Hodge3( U::RealVec , V::RealVec , W::RealVec )
    tpfl=typeof(U[1]) 

    VL = zeros(tpfl,4)
    for i=1:4
    Z  = zero(tpfl)
    for j=1:4
    for k=1:4
    for l=1:4
        Z += -ϵ(i,j,k,l)*U[j]*V[k]*W[l]
    end
    end
    end
    VL[i] = Z
    end

    return VL
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
function Hodge4( U::RealVec , V::RealVec , W::RealVec , Z::RealVec )
    tpfl=typeof(U[1]) 

    H = zero(tpfl)
    for i=1:4
    for j=1:4
    for k=1:4
    for l=1:4
        H += -ϵ(i,j,k,l)*U[i]*V[j]*W[k]*Z[l]
    end
    end
    end
    end

    return H
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
function Frame( X::RealMtx , erm::Bool=true )
    # Constructs spatial frame
    tpfl=typeof(X[1,1])       # Extracting floating point datatype
    s = 1       # Integer for constructing frame
    b = true    # Boolean variable - stops routine if not spacelike

    E = zeros(tpfl,4,3)

    for i=1:3
        E[:,i] = X[:,i] - X[:,4]
    end

    return E
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
function ConfVec( E::RealMtx )
    # Constructs configuration vector
    tpfl=typeof(E[1,1])

    NVL = Hodge3( E[:,1] , E[:,2] , E[:,3] )

    return ημν( tpfl , 4 )*NVL
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
function TrVec( E::RealMtx , thresh=1e-13 )
    tpfl=typeof(E[1,1])
    χ = ConfVec( E )

    ξ = χ

    if dot(χ,χ)>0 
        if η(ξ,χ) == 0 || abs(η(ξ,χ))<thresh
            ξ[2:4] = -χ[2:4] 
        end
        while η(ξ,χ) == 0 || abs(η(ξ,χ))<thresh
            ξ = rand(tpfl,4)
        end
        return ξ
    else
        return zeros(tpfl,4)
    end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
function ystar( E::RealMtx )
    # Constructs timelike normal to frame
    tpfl=typeof(E[1,1])

    χ = ConfVec( E )
    ξ = rand(4)

    h23 = Hodge2( E[:,2] , E[:,3] )
    h31 = Hodge2( E[:,3] , E[:,1] )
    h12 = Hodge2( E[:,1] , E[:,2] )

    Ω1 = η(E[:,1],E[:,1])/2
    Ω2 = η(E[:,2],E[:,2])/2
    Ω3 = η(E[:,3],E[:,3])/2

    H  = Ω1*h23 + Ω2*h31 + Ω3*h12

    ys = zeros(tpfl,4)

    for i=1:4
    Z  = zero(tpfl)
    for j=1:4
        Z += ξ[j]*H[j,i]
    end
    ys[i] = Z
    end

    ys = ys / η(ξ,χ)

    return ημν( tpfl , 4 )*ys
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
function locator( X::RealMtx )
    tpfl=typeof(X[1,1])

    E   = Frame( X )
    χ   = ConfVec( E )
    ys  = ystar( E )

    φ1 = η(ys,ys) / (η(ys,χ) + sqrt( abs(η(ys,χ)^2 - η(ys,ys)*η(χ,χ)) ))
    φ2 = η(ys,ys) / (η(ys,χ) - sqrt( abs(η(ys,χ)^2 - η(ys,ys)*η(χ,χ)) ))

    x1 = X[:,4] + ys - φ1*χ
    x2 = X[:,4] + ys - φ2*χ

    return (x1,x2)
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
end     # End scope of module cereal
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
module cerealtest    # Test module for cereal
#-----------------------------------------------------------------------

using LinearAlgebra
import ..cereal         # Importing functions from cereal

const RealVec{T<:Real} = Array{T,1}      # Defining vector datatype
const RealMtx{T<:Real} = Array{T,2}      # Defining matrix datatype

#-----------------------------------------------------------------------
function vgenerator( tpfl::DataType )    # Generates random 3-vector
    θ = pi*rand(tpfl)
    ϕ = 2*pi*rand(tpfl)
    V = [   sin(θ)*cos(ϕ)   ;
            sin(θ)*sin(ϕ)   ;
            cos(θ)          ]
    return V
end  # End vgenerator

#-----------------------------------------------------------------------
function nvgenerator( tpfl::DataType )   # Generates random 4-velocity
    v = vgenerator(tpfl)
    v = rand(tpfl)*rand(tpfl)*v
    V = [ one(tpfl) ; v[1] ; v[2] ; v[3] ]
    fl = one(tpfl) / sqrt(abs(cereal.η(V,V)))

    NV = fl*V
    return NV
end  # End nvgenerator

#-----------------------------------------------------------------------
function TVchk( Y::RealMtx , δV::Real=1e-14 )
    # Check tetrahedron volume is > |δV|    tpfl=typeof(Y[1,1])
    tpfl=typeof(Y[1,1])

    x = zeros(tpfl,3,4)     # Create containers
    v = zeros(tpfl,3,3)

    x[1,:] = Y[2,:]         # Spatial points for corners of tetrahedron
    x[2,:] = Y[3,:]
    x[3,:] = Y[4,:]

    v[:,1] = x[:,2] - x[:,1]    # "Frame" vectors centered on x[:,1]
    v[:,2] = x[:,3] - x[:,1]
    v[:,3] = x[:,4] - x[:,1]

    Vol = abs(det([ v[:,1]  v[:,2]  v[:,3] ]))/6    # Volume

    if Vol > abs(δV)
        return true
    else
        return false
    end
end  # End TVchk

#-----------------------------------------------------------------------
function pointgenerator( tpfl::DataType , δV::Real=1e-14 )     
    # Generates four random points W all with same time coordinate
    b = false

    W = zeros(tpfl,4,4)

    while b == false
        V1 = (one(tpfl)+rand(tpfl))*vgenerator(tpfl)/(2*one(tpfl))
        V2 = (one(tpfl)+rand(tpfl))*vgenerator(tpfl)/(2*one(tpfl))
        V3 = (one(tpfl)+rand(tpfl))*vgenerator(tpfl)/(2*one(tpfl))
        V4 = (one(tpfl)+rand(tpfl))*vgenerator(tpfl)/(2*one(tpfl))

        Vt = rand(tpfl)

        W = [   Vt     Vt     Vt     Vt     ;
                V1[1]  V2[1]  V3[1]  V4[1]  ;
                V1[2]  V2[2]  V3[2]  V4[2]  ;
                V1[3]  V2[3]  V3[3]  V4[3]  ]
        
        b = TVchk(W,δV)
    end

    return W
end  # End pointgenerator

#-----------------------------------------------------------------------
function LTM( NV::RealVec )     
    # Constructs Lorentz transformation matrix
    tpfl=typeof(NV[1])
        l = one(tpfl)           # l has a value of one
        γ = NV[1]               # γ is the Lorentz factor
        δ = (γ - l)             # A useful quantity
        NVs = NV[2]^2 + NV[3]^2 + NV[4]^2      # Square of spatial part
        if NVs > 0              # Normalization factor check
            nf = 1/sqrt(NVs)    # Normalization factor
            β = sqrt(NVs/(γ*γ)) # β is v/c
        else
            nf = zero(tpfl)     # Set normalization factor to zero
            β = zero(tpfl)      # β is v/c
        end
        vx = nf*NV[2]     # x-component of unit vector
        vy = nf*NV[3]     # y-component of unit vector
        vz = nf*NV[4]     # z-component of unit vector
    return [  γ        -γ*β*vx       -γ*β*vy       -γ*β*vz      ;
             -γ*β*vx   l + δ*(vx^2)  δ*vx*vy       δ*vx*vz      ;
             -γ*β*vy   δ*vy*vx       l + δ*(vy^2)  δ*vy*vz      ;
             -γ*β*vz   δ*vz*vx       δ*vz*vy       l + δ*(vz^2) ]
end  # End LTM

#-----------------------------------------------------------------------
function LT( W::RealMtx, NV::RealVec )     
    # Lorentz transforms W to frame where vector NV is in t direction
    tpfl=typeof(W[1,1])
    l = one(tpfl)           # l has a value of one
    Z = zeros(tpfl,4,4)

    Λ = LTM(NV)
    β = rand(tpfl)/(tpfl(2))
    γ = l/sqrt( l - β^2 )
    δ = ( γ - l )

    vx = β*NV[1]            # components of unit vector
    vy = β*NV[2]
    vz = β*NV[3]

    Λ = [  γ        -γ*β*vx        -γ*β*vy       -γ*β*vz       ;
           -γ*β*vx   l + δ*(vx^2)  δ*vx*vy       δ*vx*vz       ;
           -γ*β*vy   δ*vy*vx       l + δ*(vy^2)  δ*vy*vz       ;
           -γ*β*vz   δ*vz*vx       δ*vz*vy       l + δ*(vz^2)  ]

    for i=1:4
        Z[:,i] = Λ * W[:,i]         # Lorentz transformation
    end
    return Z
end  # End LT

#-----------------------------------------------------------------------
function epgen( tpfl::DataType=Float64 , δV::Real=1e-14 )    
    # Randomly generates emission pts
    W  = pointgenerator(tpfl,δV)
    NV = nvgenerator(tpfl)
    return LT(W,NV)
end  # End epgen

#-----------------------------------------------------------------------
function nsqr( V1::RealVec )        # Squared norm ratio
    return abs(cereal.η(V1,V1))
end  # End nsqr

#-----------------------------------------------------------------------
function single( q::Real , P::RealVec , X::RealMtx )
    # Checks if separation vectors are null
    tpfl=typeof(q)
    Y = copy(P)

    VT = [ nsqr( Y - X[:,1] ) ; nsqr( Y - X[:,2] ) ;
           nsqr( Y - X[:,3] ) ; nsqr( Y - X[:,4] ) ]

    av = ( VT[1] + VT[2] + VT[3] + VT[4] ) / (4*one(tpfl))

    if av < q
        return true
    else
        return ( X , av )
    end
end  # End single

#-----------------------------------------------------------------------
function full( iters::Number , q::Real , δV::Real=1e-14 
               , erm::Bool=true , counter::Bool=true )
    # Main test function--q determines floating point datatype
    tpfl=typeof(q)
    lb=false
    mi = zero(Int64)

    for i=1:Int64(iters)
        B = true
        if counter
            print("\r$i")
        end
        while B
            X = epgen(tpfl,δV)
            P = cereal.locator(X,erm,true)
            b = single(q,P[1],X)
            if b != true
                print("\n", b[2],"\n")
                lb = true
                mi += 1
            end
            if P[2]
                B = false
            end
        end
    end
    if lb
        print("\rTest ended with ",mi," failed cases out of "
              ,iters," \n")
    else
        print("\rTest ended with zero failed cases out of ",iters," \n")
    end
end  # End full

#-----------------------------------------------------------------------
end     # End scope of module cerealtest
#-----------------------------------------------------------------------