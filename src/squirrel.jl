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

include("metric.jl")
include("metrics/Minkowski.jl")
include("srl/FHC21.jl")
include("srl/RTC21.jl")
include("srl/mloc.jl")

#-----------------------------------------------------------------------
#
#   GEODESIC ENDPOINT AND JACOBIAN
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    gsolve( Xi::RealVec , Vi::RealVec , g::Function , tol::Real , integrator=AutoVern7(Rodas5()) )

The `gsolve` function takes an initial point `Xi` and four velocity `Vi`
and computes the endpoint of a future pointing null geodesic in the
metric func `g`. The variable `integrator` specifies the integration
scheme, and `tol` is the tolerance parameter.

"""
function gsolve( Xi::RealVec , Vi::RealVec , g::Function , tol::Real
                 , integrator=AutoVern7(Rodas5()) )
    tpfl=typeof(Xi[1])
    Z0 = zeros(tpfl,8)
    V0 = nullenforcerf( Vi , Xi , g )
    Z0 = vcat( Xi , V0 )
    return solveZ( Z0 , g , tol , tol , integrator , tol )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    V34( V3::RealVec )

The function `V34` takes a three-vector and constructs from it a
four-vector with vanishing time component.

"""
function V34( V3::RealVec )
    tpfl = typeof(V3[1])
    return [ tpfl(0) ; V3[1] ; V3[2] ; V3[3] ]
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    VidF( Zi::RealMtx )

The function `VidF` takes a matrix of initial data points and constructs
from it a single vector of initial 3-velocities for null vectors.

"""
function VidF( Zi::RealMtx )
    return vcat( Zi[6:8,1] , Zi[6:8,2] , Zi[6:8,3] , Zi[6:8,4] )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    zF( Vid::RealVec , Zi::RealMtx , g::Function , tol::Real )

The function `zF` returns differences between the endpoints of four
geodesics for the initial data encoded in `Vid` and `Zi`, and the metric
function `g`. The variable `tol` is the tolerance parameter for the
integration. This function vanishes when the four geodesics intersect.

"""
function zF( Vid::RealVec , Zi::RealMtx , g::Function , tol::Real )
    # Want to find roots of this wrt initial velocity
    tpfl=typeof(Vid[1])

    Xf = [zeros(tpfl,4) for _ in 1:4]

    Threads.@threads for i=1:4
        ( k1 , k2 ) = ( 1 + 3*(i-1) , 3 + 3*(i-1) )
        Xf[i] = gsolve( Zi[1:4,i] , V34( Vid[k1:k2] ) , g , tol )[1:4]
    end

    return vcat( Xf[1] - Xf[2] , Xf[1] - Xf[3] , Xf[1] - Xf[4] )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    gejac( Xi::RealVec , Vi::RealVec , g::Function , δ::Real )

The function `gejac` computes the endpoint of a geodesic and its
Jacobian. The variables have the same meaning as those in `gsolve`.

"""
function gejac( Xi::RealVec , Vi::RealVec , g::Function , tol::Real )
    v = Vi[2:4]
    result = DiffResults.JacobianResult(vcat(Xi,Vi),v)
    result = ForwardDiff.jacobian!(result,var->gsolve(Xi,V34(var),g,tol)
                                    , v )
    return ( DiffResults.value(result)
             , DiffResults.jacobian(result)[1:4,:] )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    geocJ( Zi::RealMtx , g::Function , δ::Real )

The function `geocJ` computes the Jacobian of the function `zF` from the
Jacobian endpoints computed in `gejac`.

"""
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
"""
    zFc( Zf::RealMtx )

The function `zFc` rearranges the endpoints `Zf` to match the output of
the function `zF`.

"""
function zFc( Zf::RealMtx )
    return vcat( Zf[1:4,1] - Zf[1:4,2] , Zf[1:4,1] - Zf[1:4,3] ,
                 Zf[1:4,1] - Zf[1:4,4] )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    idf( Zi::RealMtx , gfunc::Function , tol::Real , nb::Int )

The function `idf` takes as input `Zi`, a matrix formed from the
emission points and guesses for the initial four-velocities, and the
metric function `gfunc`. It outputs the corrected initial data for
intersecting null geodesics. The argument `tol` is the tolerance
parameter and `nb` is the Broyden termination limit. This function calls
`geocJ` to compute the initial Jacobian, then calls and applies the
Broyden root finding function `bsolve` to the function `zF` (also
called).

"""
function idf( Zi::RealMtx , gfunc::Function , tol::Real , nb::Int )
    tpfl=typeof(Zi[1,1])

    Zf = copy(Zi)
    V0 = VidF(Zi)
    res = geocJ( Zi , gfunc , tol )

    ( F0 , J ) = ( zFc(res[1]) , res[2] )

    V = bsolve( v->zF(v,Zi,gfunc,tol) , J , F0 , V0 , nb )

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
"""
    locator4( X::RealMtx , Xc::RealVec , gfunc::Function , tol::Real , nb::Int=24 , idv::Bool=false , V::RealMtx=zeros(Float64,4,4) )

The function `locator4` computes the intersection point from a set of
four emission points `X` using the guess `Xc`. The intersection point is
computed by first finding the initial data using the function `idf`,
then integrating geodesics to obtain the result. If an improved guess
for the four-velocities is available, one can set `ìdv=true` and specify
the four-velocities as column vectors in the matrix `V`.

"""
function locator4( X::RealMtx , Xc::RealVec , gfunc::Function ,
                   tol::Real , nb::Int=24 , idv::Bool=false , 
                   V::RealMtx=zeros(Float64,4,4) )
    tpfl=typeof(X[1,1])

    Zi  = zeros(tpfl,8,4)
    Zf = copy(Zi)

    if idv 
        for i=1:4
            Zi[:,i] = vcat( X[:,i] , V[:,i] )
        end
    else
        for i=1:4
            Zi[:,i] = vcat( X[:,i] , Xc - X[:,i] )
        end
    end

    Zi = idf( Zi , gfunc , tol , nb )
    Threads.@threads for i=1:4
        Zf[:,i] = gsolve( Zi[1:4,i] , Zi[5:8,i] , gfunc , tol )
    end
    return (Zf[1:4,1] + Zf[1:4,2] + Zf[1:4,3] + Zf[1:4,4])/4 
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#
#   RELATIVISTIC LOCATION WITH >4 EMISSION POINTS
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    locator( X::RealMtx , gfunc::Function , δ::Real , ne::Int=5 , nb::Int=24 , outthresh::Real=2e1 , tpflc::DataType=Double64 )

The function `locator` computes the intersection point from a set of
`ne>4` emission points `X` by applying `locator4` to all combinations of
4 points out of `ne` in `X`. A basic outlier detection algorithm
(implemented in the function `odetc`) is applied to reduce errors.

"""
function locator(  X::RealMtx , gfunc::Function , δ::Real , ne::Int=5 ,
                   nb::Int=24 , outthresh::Real=2e1 , 
                   tpflc::DataType=Double64 )
    tpfl  = typeof(X[1,1])

    l = size(X)

    if gfunc == x->ημν(x)
        return mlocator( tpflc.(X) )
    elseif  l[2] < 4 || ne < 4
        print("Need more than four emission points.")
        return zeros(tpfl,4)
    elseif l[2] == 4 || ne == 4 
        Xdual   = locator4FHC21( tpflc.(X) )
        X1 = locator4( X , Xdual[1] , gfunc , δ , nb , false )
        X2 = locator4( X , Xdual[2] , gfunc , δ , nb , false )
        return (X1,X2)
    elseif l[2] >= 5
        if  ne >= 5 && ne < size(X)[2]
            Y = X[:,1:ne]
        elseif size(X)[2] == 5 || ne >= size(X)[2]
            Y = X
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
