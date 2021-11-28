#-----------------------------------------------------------------------
module squirrel     # BEGIN MODULE SQUIRREL
#-----------------------------------------------------------------------

using LinearAlgebra, Combinatorics, Statistics , DoubleFloats
using ForwardDiff, DiffResults

include("cereal.jl")
include("type.jl")
include("geosol.jl")
include("broyden.jl")
include("outlier.jl")

#-----------------------------------------------------------------------
#
#   GEODESIC INTEGRATOR AND ZERO FUNCTION
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#   GEODESIC SOLVER
#-----------------------------------------------------------------------
function gsolve( Xi::RealVec , Vi::RealVec , gfunc::Function , δ::Real 
                 , integrator=AutoVern7(Rodas5()) )
    tpfl=typeof(Xi[1])
    Z0 = zeros(tpfl,8)
    V0 = nullenforcerf( Vi , Xi , gfunc )
    Z0 = vcat( Xi , V0 )
    return solveZ( Z0 , gfunc , δ , δ , integrator , δ )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   ADD COMPONENT OF VECTOR
#-----------------------------------------------------------------------
function V34( V3::RealVec )
    tpfl = typeof(V3[1])
    return [ tpfl(0) ; V3[1] ; V3[2] ; V3[3] ]
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   VELOCITY INITIAL DATA FUNCTION
#-----------------------------------------------------------------------
function VidF( Zi::RealMtx )
    return vcat( Zi[6:8,1] , Zi[6:8,2] , Zi[6:8,3] , Zi[6:8,4] )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   ZERO FUNCTION
#-----------------------------------------------------------------------
function zF( Vid::RealVec , Zi::RealMtx , gfunc::Function , δ::Real )
    # Want to find roots of this wrt initial velocity
    tpfl=typeof(Vid[1])

    Xf = [zeros(tpfl,4) for _ in 1:4]

    Threads.@threads for i=1:4
        ( k1 , k2 ) = ( 1 + 3*(i-1) , 3 + 3*(i-1) )
        Xf[i] = gsolve( Zi[1:4,i] , V34( Vid[k1:k2] ) , gfunc , δ )[1:4]
    end

    return vcat( Xf[1] - Xf[2] , Xf[1] - Xf[3] , Xf[1] - Xf[4] )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#
#   GEODESIC ENDPOINT AND JACOBIAN
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#   CALCULATE GEODESIC ENDPOINT AND ITS JACOBIAN
#-----------------------------------------------------------------------
function gejac( X::RealVec , V::RealVec , g::Function , δ::Real )
    v = V[2:4]
    result = DiffResults.JacobianResult(vcat(X,V),v)
    result = ForwardDiff.jacobian!( result , var->gsolve(X,V34(var),g,δ)
                                    , v )
    return ( DiffResults.value(result)
             , DiffResults.jacobian(result)[1:4,:] )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   JACOBIAN CALCULATOR
#-----------------------------------------------------------------------
function geocJ( Zi::RealMtx , g::Function , δ::Real )
    tpfl=typeof(Zi[1,1])
    v   = zeros(tpfl,3)
    J   = zeros(tpfl,12,12)
    dXV = [zeros(tpfl,4,3) for _ in 1:4]
    Zf  = copy(Zi)

    Threads.@threads for i=1:4
        res = gejac( Zi[1:4,i] , Zi[5:8,i] , g , δ )
        ( Zf[:,i] , dXV[i] ) = ( res[1] , res[2] )
    end

    J[1:4,1:3]  = dXV[1]
    J[5:8,1:3]  = dXV[1]
    J[9:12,1:3] = dXV[1]

    J[1:4,4:6]    = - dXV[2]
    J[5:8,7:9]    = - dXV[3]
    J[9:12,10:12] = - dXV[4]

    return (Zf,J)
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#
#   SOLVING FOR INITIAL DATA
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#   REARRANGE ENDPOINTS
#-----------------------------------------------------------------------
function zFc( Zf::RealMtx )
    return vcat( Zf[1:4,1] - Zf[1:4,2] , Zf[1:4,1] - Zf[1:4,3] ,
                 Zf[1:4,1] - Zf[1:4,4] )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   INITIAL DATA CORRECTOR
#-----------------------------------------------------------------------
function idc( Zi::RealMtx , gfunc::Function , δ::Real , nb::Int )
    tpfl=typeof(Zi[1,1])

    Zf = copy(Zi)
    V0 = VidF(Zi)
    res = geocJ( Zi , gfunc , δ )

    ( F0 , J ) = ( zFc(res[1]) , res[2] )

    V = bsolve( v->zF(v,Zi,gfunc,δ) , J , F0 , V0 , nb )

    Zf[5:8,1] = V34(  V[1:3]  )
    Zf[5:8,2] = V34(  V[4:6]  )
    Zf[5:8,3] = V34(  V[7:9]  )
    Zf[5:8,4] = V34( V[10:12] )

    return Zf
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#
#   RELATIVISTIC LOCATOR
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#   SINGLE LOCATOR
#-----------------------------------------------------------------------
function slocator( X::RealMtx , gfunc::Function , δ::Real , nb::Int=24 
                   , alt::Bool=false , Xg::RealVec=zeros(Float64,4) 
                   , erc::Bool=false )
    tpfl=typeof(X[1,1])

    Zi  = zeros(tpfl,8,4)

    if alt
        Zf = copy(Zi)
        Xc = Xg

        for i=1:4
            Zi[:,i] = vcat( X[:,i] , Xc - X[:,i] )
        end
        Zi = idc( Zi , gfunc , δ , nb )
        Threads.@threads for i=1:4
            Zf[:,i] = gsolve( Zi[1:4,i] , Zi[5:8,i] , gfunc , δ )
        end
        return ( (Zf[1:4,1] + Zf[1:4,2] + Zf[1:4,3] + Zf[1:4,4])/4 ,
                 zeros(tpfl,4) )
    else
        ZfA = copy(Zi)
        ZfB = copy(Zi)

        Xc = cereal.slocator(X,false)

        for i=1:4
            Zi[:,i] = vcat( X[:,i] , Xc[1] - X[:,i] )
        end
        Zi = idc( Zi , gfunc , δ , nb )
        Threads.@threads for i=1:4
            ZfA[:,i] = gsolve( Zi[1:4,i] , Zi[5:8,i] , gfunc , δ )
        end
    
        for i=1:4
            Zi[:,i] = vcat( X[:,i] , Xc[2] - X[:,i] )
        end
        Zi = idc( Zi , gfunc , δ , nb )
        Threads.@threads for i=1:4
            ZfB[:,i] = gsolve( Zi[1:4,i] , Zi[5:8,i] , gfunc , δ )
        end
    
        return ( (ZfA[1:4,1] + ZfA[1:4,2] + ZfA[1:4,3] + ZfA[1:4,4])/4 ,
                 (ZfB[1:4,1] + ZfB[1:4,2] + ZfB[1:4,3] + ZfB[1:4,4])/4 )
    end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#
#   RELATIVISTIC LOCATION WITH >4 EMISSION POINTS
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#   MULTI LOCATOR
#-----------------------------------------------------------------------
function mlocator( Yi::RealMtx , gfunc::Function , δ::Real ,
                   nb::Int=24 , ξ1::Real=Double64(1e-18) , 
                   ξ2::Real=1e1 , ne::Int=6 )
    tpfl=typeof(Yi[1,1])

    if  size(Yi)[2] <= 4 && ne <= 4
        print("Need more than four emission points.")
        return zeros(tpfl,4)
    else
        if  ne > 4 && ne < size(Yi)[2]
            Y = Yi[:,1:ne]
        elseif ne >= size(Yi)[2]
            Y = Yi
        end

        Xc      = cereal.mlocator( Y , Double64(ξ1) , false )[1]
        W       = combX(Y)
        nr      = length(W)
        Xs      = [zeros(tpfl,4) for _ in 1:nr]
        Xsa     = [zeros(tpfl,4) for _ in 1:nr]

        for i=1:nr
            Xs[i]  = slocator( W[i] , gfunc , δ , nb , true ,
                               tpfl.(Xc) , false )[1]
        end

        return mean(odetc( Xs , ξ2 )[1])
    end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
end     # END MODULE SQUIRREL
#-----------------------------------------------------------------------
