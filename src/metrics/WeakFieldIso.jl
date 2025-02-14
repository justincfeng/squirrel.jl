#-----------------------------------------------------------------------
#   POTENTIAL FUNCTIONS
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#   WGS84 Earth reference potential function
#-----------------------------------------------------------------------
"""
	Vpot( X::RealVec , p::RealVec )
The `Vpot` function takes a point `X` in Cartesian coordinates and 
calculates the gravitational potential of an Earth-like object with parameters:
- p[1]: GM (product of gravitational constant and mass)
- p[2]: J2 (quadrupole moment)
- p[3]: a (equatorial radius)
"""
function Vpot( X::RealVec , p::RealVec )
    tpfl=typeof(X[1])
    GM = p[1]
    J2 = p[2]
    a = p[3]
    rs = dot(X[2:4],X[2:4])
    return -( GM/sqrt(rs) )*( tpfl(1) - J2*((a)^2/rs)*Pl(cθ(X),2) )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   WGS84 Earth reference ellipsoid radius
#-----------------------------------------------------------------------
"""
	rell( x::RealVec , p::RealVec )
The `rell` function takes a direction defined by vector `x` and 
calculates the distance along that direction from the origin to the 
surface of an oblate spheroidal ellipsoid with parameters:
- p[1]: a (semimajor axis)
- p[2]: b (semiminor axis)
"""
function rell( x::RealVec , p::RealVec )
    tpfl=typeof(x[1])
    a = p[1]
    b = p[2]
    rs = dot(x,x)
    return tpfl(a*b)/sqrt( b^2 + (a - b)*(a + b)*(x[3]^2)/rs )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   WGS84 Earth reference ellipsoid rescale factor
#-----------------------------------------------------------------------
"""
	rellsc( x::RealVec , p::RealVec )
The `rellsc` function takes a direction defined by vector `x` and 
calculates the relative distance (divided by a) along that direction 
from the origin to the surface of an oblate spheroidal ellipsoid with parameters:
- p[1]: a (semimajor axis)
- p[2]: b (semiminor axis)
"""
function rellsc( x::RealVec , p::RealVec )
    tpfl=typeof(x[1])
    a = p[1]
    b = p[2]
    rs = dot(x,x)
    return tpfl(b)/sqrt( b^2 + (a - b)*(a + b)*(x[3]^2)/rs )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   ISOTROPIC METRIC COMPONENTS
#-----------------------------------------------------------------------
"""
	giso( X::RealVec , p::RealVec )
The `giso` function takes a point `X` in Cartesian coordinates and 
calculates the weak-field metric components for an Earth-like object with parameters:
- p[1]: GM (product of gravitational constant and mass)
- p[2]: J2 (quadrupole moment)
- p[3]: a (equatorial radius)
"""
function giso( X::RealVec , p::RealVec )
    tpfl=typeof(X[1])
    return ημν(tpfl) + Vpot(X,p)*δμν(tpfl)
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    ge_iso( X::RealVec )
The `ge_iso` function takes a point `X` in Cartesian coordinates and 
calculates the weak-field metric components for Earth with default parameters:
- GM = 1 (product of gravitational constant and mass)
- J2 = 1.0826300e-3 (quadrupole moment)
- a = 1.438127773656399e9 (equatorial radius)
"""
function ge_iso( X::RealVec )
    tpfl=typeof(X[1])
    p = [tpfl(1), tpfl(1.0826300e-3), tpfl(1.438127773656399e9)]  # [GM, J2, a]
    return giso(X,p)
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
