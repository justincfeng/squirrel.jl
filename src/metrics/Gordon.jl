#-----------------------------------------------------------------------
#       BEGIN   Gordon.jl
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#   DYAD
#-----------------------------------------------------------------------
"""
    uugen( u::RealVec, X::RealVec )

The function `uugen` constructs a dyad from the vector u.

"""
function uugen( u::RealVec, X::RealVec )
    return transpose(u).*u
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   COSINE OF THETA
#-----------------------------------------------------------------------
"""
    cθ( X::RealVec )

The function `cθ` computes the cosine of θ (as defined in spherical
polar coordinates, using the physicist convention).

"""
function cθ( X::RealVec )
    return X[4]/norm(X[2:4])
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   SINE OF PHI
#-----------------------------------------------------------------------
"""
    sϕ( X::RealVec )

The function `sϕ` computes the sine of ϕ (as defined in spherical
polar coordinates, using the physicist convention).

"""
function sϕ( X::RealVec )
    return X[3]/norm(X[3:4])
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   THE GORDON METRIC
#-----------------------------------------------------------------------
"""
    gGordon( X::RealVec , p::RealVec , n::Function=nIR , gfunc::Function=giso , 
             U::RealVec=Float64[-1;0;0;0] )

The function `gGordon` computes the Gordon metric. The arguments are:
- X: coordinates
- p: parameter vector passed to the background metric function gfunc
- n: index of refraction function (can be replaced with user-supplied function)
- gfunc: background metric function (can be replaced with user-supplied function)
- U: fluid four-velocity function (can be replaced with user-supplied vector)
"""
function gGordon( X::RealVec , p::RealVec , n::Function=nIR ,
                  gfunc::Function=giso , U::RealVec=Float64[-1;0;0;0] )

    tpfl=typeof(X[1])

    gs = gfunc( X, p )

    Unsq = abs(transpose(U)*gs*U)

    UU = uugen( U , X )

    return gs + ( 1-1/(n(X)^2) )*UU/Unsq
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    ge_gordon( X::RealVec )
The `ge_gordon` function computes the Gordon metric with default Earth parameters:
- GM = 1 (product of gravitational constant and mass)
- J2 = 1.0826300e-3 (quadrupole moment)
- a = 1.438127773656399e9 (equatorial radius)
"""
function ge_gordon( X::RealVec )
    tpfl=typeof(X[1])
    p = [tpfl(1), tpfl(1.0826300e-3), tpfl(1.438127773656399e9)]  # [GM, J2, a]
    return gGordon(X, p)
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#       END     Gordon.jl
#-----------------------------------------------------------------------