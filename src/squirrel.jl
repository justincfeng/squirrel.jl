#-----------------------------------------------------------------------
module squirrel     # BEGIN MODULE SQUIRREL
#-----------------------------------------------------------------------

using LinearAlgebra, Combinatorics, Statistics , DoubleFloats
using ForwardDiff, DiffResults

include("type.jl")
include("seval.jl")
include("geosol.jl")
include("broyden.jl")
include("outlier.jl")

include("metrics/Minkowski.jl")
include("srl/FHC21.jl")
include("srl/RTC21.jl")
include("srl/mloc.jl")

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
#   LOCATOR FOR FOUR EMISSION POINTS
#-----------------------------------------------------------------------
function locator4( X::RealMtx , Xc::RealVec , gfunc::Function ,
                   δ::Real , nb::Int=24 , erc::Bool=false )
    tpfl=typeof(X[1,1])

    Zi  = zeros(tpfl,8,4)
    Zf = copy(Zi)

    for i=1:4
        Zi[:,i] = vcat( X[:,i] , Xc - X[:,i] )
    end
    Zi = idc( Zi , gfunc , δ , nb )
    Threads.@threads for i=1:4
        Zf[:,i] = gsolve( Zi[1:4,i] , Zi[5:8,i] , gfunc , δ )
    end
    return (Zf[1:4,1] + Zf[1:4,2] + Zf[1:4,3] + Zf[1:4,4])/4 
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#
#   RELATIVISTIC LOCATION WITH >4 EMISSION POINTS
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#   MULTI LOCATOR
#-----------------------------------------------------------------------
function locator(  Yi::RealMtx , gfunc::Function , δ::Real ,
                   nb::Int=24 , tpflc::DataType=Double64 , 
                   outthresh::Real=1e1 , ne::Int=5 )
    tpfl  = typeof(Yi[1,1])

    l = size(Yi)

    if  l[2] < 4 || ne < 4
        print("Need more than four emission points.")
        return zeros(tpfl,4)
    elseif l[2] == 4 || ne == 4 
        Xdual   = locator4FHC21( tpflc.(Y) )
        X1 = locator4( Yi , Xdual[1] , gfunc , δ , nb , false )
        X2 = locator4( Yi , Xdual[1] , gfunc , δ , nb , false )
        return (X1,X2)
    elseif l[2] >= 5
        if  ne >= 5 && ne < size(Yi)[2]
            Y = Yi[:,1:ne]
        elseif size(Yi)[2] == 5 || ne >= size(Yi)[2]
            Y = Yi
        end

        Xc      = mlocator( tpflc.(Y) )
        W       = combX(Y)
        nr      = length(W)
        Xs      = [zeros(tpfl,4) for _ in 1:nr]

        for i=1:nr
            Xs[i]  = locator4( W[i] , tpfl.(Xc) , gfunc , δ , nb , 
                               false )
        end

        return mean(odetc( Xs , outthresh )[1])
    end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
end     # END MODULE SQUIRREL
#-----------------------------------------------------------------------
