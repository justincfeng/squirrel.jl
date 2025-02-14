#-----------------------------------------------------------------------
#       BEGIN   KerrSchild.jl
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#   KERR-SCHILD FUNCTIONS
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    rsq( X::RealVec , p::RealVec )
The `rsq` function takes a point `X` in Cartesian Kerr-Schild 
coordinates and calculates the value of ``r^2`` at that point in a
Kerr spacetime with parameters `p` where p[1] is the rotation parameter a
"""
function rsq( X::RealVec , p::RealVec )
    tpfl=typeof(X[1])
    x = X[2]
    y = X[3]
    z = X[4]
    a = p[1]  # rotation parameter

    return ( -a^2 + x^2 + y^2 + z^2 + sqrt(4*(a^2)*(z^2) 
             + (-a^2 + x^2 + y^2 + z^2)^2))/2
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    fks( X::RealVec , p::RealVec )
The `fks` function takes a point `X` in Cartesian Kerr-Schild 
coordinates and calculates the value of function ``f`` at that point in 
a Kerr spacetime with parameters `p` where p[1] is the rotation parameter a
and p[2] is GM (product of gravitational constant and mass)
"""
function fks( X::RealVec , p::RealVec )
    tpfl=typeof(X[1])
    x = X[2]
    y = X[3]
    z = X[4]
    a = p[1]   # rotation parameter
    GM = p[2]  # mass parameter

    rs = rsq(X,p)
    r = sqrt(rs)
    return 2*GM*rs*r/(rs^2+(a^2)*(z^2))
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    kks( X::RealVec , p::RealVec )
The `kks` function takes a point `X` in Cartesian Kerr-Schild 
coordinates and calculates the tensor product ``k_\\mu k_\\nu`` at that 
point in a Kerr spacetime with parameters `p` where p[1] is the rotation parameter a
"""
function kks( X::RealVec , p::RealVec )
    tpfl=typeof(X[1])
    k = zeros(tpfl,4)
    kk = zeros(tpfl,4,4)
    x = X[2]
    y = X[3]
    z = X[4]
    a = p[1]  # rotation parameter

    rs = rsq(X,p)
    r = sqrt(rs)

    k[1] = one(tpfl)
    k[2] = (r*x+a*y)/(rs+a^2)
    k[3] = (r*y-a*x)/(rs+a^2)
    k[4] = z/r

    for i=1:4
        for j=1:4
            kk[i,j] = k[i]*k[j]
        end
    end
    return kk
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   KERR-SCHILD METRIC COMPONENTS
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    gks( X::RealVec , p::RealVec )
The `gks` function takes a point `X` in Cartesian Kerr-Schild 
coordinates and calculates the components of the Kerr-Schild metric at 
that point in a Kerr spacetime with parameters `p` where:
- p[1] is the rotation parameter a
- p[2] is GM (product of gravitational constant and mass)
"""
function gks( X::RealVec , p::RealVec )
    tpfl=typeof(X[1])
    return ημν(tpfl) + fks(X,p)*kks(X,p)
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    ge( X::RealVec )
The `ge` function takes a point `X` in Cartesian Kerr-Schild 
coordinates and calculates the components of the Kerr-Schild metric at 
that point in a Kerr spacetime with default Earth parameters:
- rotation parameter a = 738 (roughly corresponding to Earth's angular momentum)
- GM = 1 (product of gravitational constant and mass)
"""
function ge( X::RealVec )
    tpfl=typeof(X[1])
    p = [tpfl(738), tpfl(1)]  # [a, GM]
    return gks(X,p)
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#       END     KerrSchild.jl
#-----------------------------------------------------------------------
