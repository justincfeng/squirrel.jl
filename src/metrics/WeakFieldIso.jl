#-----------------------------------------------------------------------
#   POTENTIAL FUNCTIONS
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#   WGS84 Earth reference potential function
#-----------------------------------------------------------------------
"""
	Vpot( X::RealVec , GM::Real=1 , J2::Real=1.0826300e-3 
	, a::Real=1.438127773656399e9 )
The `Vpot` function takes a point `X` in Cartesian coordinates and 
calculates the gravitational potential of an Earth-like object with
product of gravitational constant and mass `GM`, quadrupole moment `J2`
and equatorial radius `a`
"""
function Vpot( X::RealVec , GM::Real=1 , J2::Real=1.0826300e-3 
    , a::Real=1.438127773656399e9 )
tpfl=typeof(X[1])
rs = dot(X[2:4],X[2:4])
return -( GM/sqrt(rs) )*( tpfl(1) - J2*((a)^2/rs)*Pl(cθ(X),2) )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   WGS84 Earth reference ellipsoid radius
#-----------------------------------------------------------------------
"""
	rell( x::RealVec , a::Real=1.438127773656399e9 
	, b::Real=1.433306003519573e9 )
The `rell` function takes a direction defined by vector `x` and 
calculates the distance along that direction from the origin to the 
surface of an oblate spheroidal ellipsoid with semimajor axis `a` and 
semiminor axis `b`
"""
function rell( x::RealVec , a::Real=1.438127773656399e9 
, b::Real=1.433306003519573e9 )
tpfl=typeof(x[1])
rs = dot(x,x)
return tpfl(a*b)/sqrt( b^2 + (a - b)*(a + b)*(x[3]^2)/rs )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   WGS84 Earth reference ellipsoid rescale factor
#-----------------------------------------------------------------------
"""
	rellsc( x::RealVec , a::Real=1.438127773656399e9 
	, b::Real=1.433306003519573e9 )
The `rellsc` function takes a direction defined by vector `x` and 
calculates the relative distance (divided by `a`) along that direction 
from the origin to the surface of an oblate spheroidal ellipsoid with 
semimajor axis `a` and semiminor axis `b`
"""
function rellsc( x::RealVec , a::Real=1.438127773656399e9 
, b::Real=1.433306003519573e9 )
tpfl=typeof(x[1])
rs = dot(x,x)
return tpfl(b)/sqrt( b^2 + (a - b)*(a + b)*(x[3]^2)/rs )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   ISOTROPIC METRIC COMPONENTS
#-----------------------------------------------------------------------
"""
	giso( X::RealVec )
The `giso` function takes a point `X` in Cartesian coordinates and 
calculates the the weak-field metric components for an Earth-like object
with product of gravitational constant and mass ``GM=1``, quadrupole 
moment ``J2=1.0826300e-3\times10^{-3}`` and equatorial radius 
``a_{\rm ell}=1.438127773656399\times10^9``
"""
function giso( X::RealVec )
tpfl=typeof(X[1])

V = Vpot( X, tpfl(1) , tpfl(1.0826300e-3) 
, tpfl(1.438127773656399e9) )

return ημν(tpfl) + V*δμν(tpfl)
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
