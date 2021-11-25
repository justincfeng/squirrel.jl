#-----------------------------------------------------------------------------------------------------------------------------------------------
module geocoord     # BEGIN MODULE GEOCOORD
#-----------------------------------------------------------------------------------------------------------------------------------------------

using LinearAlgebra
using DoubleFloats
using Geodesy
using CoordinateTransformations

include("type.jl")

#-----------------------------------------------------------------------------------------------------------------------------------------------
#   TRANSFORM CARTESIAN TO ECI
#-----------------------------------------------------------------------------------------------------------------------------------------------
function Tc2s( X::RealVec )
    tpfl=typeof(X[1])
    xs = SphericalFromCartesian()(X[2:4])
    return tpfl[ X[1] ; xs.r ; xs.θ ; xs.ϕ ]
end

#-----------------------------------------------------------------------------------------------------------------------------------------------
#   TRANSFORM ECI TO CARTESIAN
#-----------------------------------------------------------------------------------------------------------------------------------------------
function Ts2c( X::RealVec )
    tpfl=typeof(X[1])
    xc = CartesianFromSpherical()(Spherical(X[2],X[3],X[4]))
    return tpfl[ X[1] ; xc[1] ; xc[2] ; xc[3] ]
end

#-----------------------------------------------------------------------------------------------------------------------------------------------
#   TRANSFORM X MATRIX CARTESIAN TO ECI
#-----------------------------------------------------------------------------------------------------------------------------------------------
function T4Xc2s( Xc::RealMtx )
    tpfl=typeof(Xc[1,1])

    Xs = copy(Xc)

    Xs[:,1] = Tc2s( Xc[:,1] )
    Xs[:,2] = Tc2s( Xc[:,2] )
    Xs[:,3] = Tc2s( Xc[:,3] )
    Xs[:,4] = Tc2s( Xc[:,4] )

    return Xs
end

#-----------------------------------------------------------------------------------------------------------------------------------------------
#   TRANSFORM X MATRIX ECI TO CARTESIAN
#-----------------------------------------------------------------------------------------------------------------------------------------------
function T4Xs2c( Xs::RealMtx )
    tpfl=typeof(Xs[1,1])

    Xc = copy(Xs)

    Xc[:,1] = Ts2c( Xs[:,1] )
    Xc[:,2] = Ts2c( Xs[:,2] )
    Xc[:,3] = Ts2c( Xs[:,3] )
    Xc[:,4] = Ts2c( Xs[:,4] )

    return Xc
end

#-----------------------------------------------------------------------------------------------------------------------------------------------
#   TRANSFORM GENERATED DATA CARTESIAN TO ECI
#-----------------------------------------------------------------------------------------------------------------------------------------------
function tgenc2s( X::RealMtx , Xcx::RealVec )
    return ( T4Xc2s( X ) , Tc2s( Xcx ) )
end

#-----------------------------------------------------------------------------------------------------------------------------------------------
#   TRANSFORM GENERATED DATA ECI TO CARTESIAN
#-----------------------------------------------------------------------------------------------------------------------------------------------
function tgens2c( X::RealMtx , Xcx::RealVec )
    return ( T4Xs2c( X ) , Ts2c( Xcx ) )
end

#-----------------------------------------------------------------------------------------------------------------------------------------------
#   
#-----------------------------------------------------------------------------------------------------------------------------------------------
function Cart2LLA( X::RealVec )
    tpfl=typeof(X[1])

    um   = tpfl(4.435028039117671e-3) ;

    t = zeros(tpfl,4)

    xlla = LLAfromECEF(wgs84)(ECEF(X[2]*um,X[3]*um,X[4]*um))

    t = X[1]
	
    return (xlla,t)
end

function FourCart2LLA( X::RealMtx )
    tpfl=typeof(X[1,1])

    um   = tpfl(4.435028039117671e-3) ;

    t = zeros(tpfl,4)

    xlla1 = LLAfromECEF(wgs84)(ECEF(X[2,1]*um,X[3,1]*um,X[4,1]*um))
    xlla2 = LLAfromECEF(wgs84)(ECEF(X[2,2]*um,X[3,2]*um,X[4,2]*um))
    xlla3 = LLAfromECEF(wgs84)(ECEF(X[2,3]*um,X[3,3]*um,X[4,3]*um))
    xlla4 = LLAfromECEF(wgs84)(ECEF(X[2,4]*um,X[3,4]*um,X[4,4]*um))

    ( t[1] , t[2] , t[3] , t[4] ) = ( X[1,1] , X[1,2] , X[1,3] , X[1,4] )
	
    return (xlla1,xlla2,xlla3,xlla4,t)
end

function LLA2X( xl1 , xl2 , xl3 , xl4 , t::RealVec )
    tpfl=typeof(t[1])

    X = zeros(tpfl,4,4)

    nm   = tpfl(1/4.435028039117671e-3) ;

    xe1 = ECEFfromLLA(wgs84)(xl1)
    xe2 = ECEFfromLLA(wgs84)(xl2)
    xe3 = ECEFfromLLA(wgs84)(xl3)
    xe4 = ECEFfromLLA(wgs84)(xl4)

    ( X[1,1] , X[1,2] , X[1,3] , X[1,4] ) = ( t[1] , t[2] , t[3] , t[4] )

    X[2:4,1] = tpfl[ xe1[1]*nm ; xe1[2]*nm ;xe1[3]*nm ]
    X[2:4,2] = tpfl[ xe2[1]*nm ; xe2[2]*nm ;xe2[3]*nm ]
    X[2:4,3] = tpfl[ xe3[1]*nm ; xe3[2]*nm ;xe3[3]*nm ]
    X[2:4,4] = tpfl[ xe4[1]*nm ; xe4[2]*nm ;xe4[3]*nm ]

	return X
end

#-----------------------------------------------------------------------------------------------------------------------------------------------
end     # END MODULE GEOCOORD
#-----------------------------------------------------------------------------------------------------------------------------------------------