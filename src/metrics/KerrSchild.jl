#-----------------------------------------------------------------------
#       BEGIN   KerrSchild.jl
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#   KERR-SCHILD FUNCTIONS
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    rsq( X::RealVec , a::Real )
The `rsq` function takes a point `X` in Cartesian Kerr-Schild 
coordinates and calculates the value of ``r^2`` at that point in a
Kerr spacetime with rotation parameter `a`
"""
function rsq( X::RealVec , a::Real )
    tpfl=typeof(X[1])
    x = X[2]
    y = X[3]
    z = X[4]

    return ( -a^2 + x^2 + y^2 + z^2 + sqrt(4*(a^2)*(z^2) 
             + (-a^2 + x^2 + y^2 + z^2)^2))/2
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    fks( X::RealVec , a::Real , GM::Real=1 )
The `fks` function takes a point `X` in Cartesian Kerr-Schild 
coordinates and calculates the value of function ``f`` at that point in 
a Kerr spacetime with rotation parameter `a` and with the product of
gravitational constant and mass `GM`
"""
function fks( X::RealVec , a::Real , GM::Real=1 )
    tpfl=typeof(X[1])
    x = X[2]
    y = X[3]
    z = X[4]

    rs = rsq(X,a)
    r = sqrt(rs)
    return 2*GM*rs*r/(rs^2+(a^2)*(z^2))
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    kks( X::RealVec , a::Real )
The `kks` function takes a point `X` in Cartesian Kerr-Schild 
coordinates and calculates the tensor product ``k_\\mu k_\\nu`` at that 
point in a Kerr spacetime with rotation parameter `a`
"""
function kks( X::RealVec , a::Real )
    tpfl=typeof(X[1])
    k = zeros(tpfl,4)
    kk = zeros(tpfl,4,4)
    x = X[2]
    y = X[3]
    z = X[4]

    rs = rsq(X,a)
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
    gks( X::RealVec , a::Real=0 , GM::Real=1 )
The `gks` function takes a point `X` in Cartesian Kerr-Schild 
coordinates and calculates the components of the Kerr-Schild metric at 
that point in a Kerr spacetime with rotation parameter `a` and with the 
product of gravitational constant and mass `GM`
"""
function gks( X::RealVec , a::Real=0 , GM::Real=1 )
    tpfl=typeof(X[1])

    return ημν(tpfl) + fks(X,a,GM)*kks(X,a)
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    ge( X::RealVec )
The `ge` function takes a point `X` in Cartesian Kerr-Schild 
coordinates and calculates the components of the Kerr-Schild metric at 
that point in a Kerr spacetime with relation between mass and angular 
momentum roughly corresponding to the Earth value (``a=738GM``)
"""
function ge( X::RealVec )
    tpfl=typeof(X[1])

    a = tpfl(738)   #   This is the roughly the rotation parameter 
                    #   for Earth in units where M_Earth = 1
                    #   For reference, the avg radius of Earth is
                    #   1.437e9

    return gks(X,a,1)  
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#       END     KerrSchild.jl
#-----------------------------------------------------------------------
