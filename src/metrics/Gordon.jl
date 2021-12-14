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
    gGordon( X::RealVec , n::Function=nIR , gfunc::Function=giso , 
             U::RealVec=Float64[-1;0;0;0] )

The function `gGordon` computes the Gordon metric. The argument `n` can
be replaced with a user-supplied index of refraction function, the
argument `gfunc` can be replaced with a user-supplied background metric,
and the argument `U` can be replaced with a user-supplied fluid 
four-velocity function.

"""
function gGordon( X::RealVec , n::Function=nIR ,
                  gfunc::Function=giso , U::RealVec=Float64[-1;0;0;0] )

    tpfl=typeof(X[1])

    gs = gfunc( X )

    Unsq = abs(transpose(U)*gs*U)

    UU = uugen( U , X )

    return gs + ( 1-1/(n(X)^2) )*UU/Unsq
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#       END     Gordon.jl
#-----------------------------------------------------------------------