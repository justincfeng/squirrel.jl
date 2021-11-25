#-----------------------------------------------------------------------
#       BEGIN   Gordon.jl
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#   DYAD
#-----------------------------------------------------------------------
function uugen( u::RealVec, X::RealVec )
    return transpose(u).*u
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   COSINE OF THETA
#-----------------------------------------------------------------------
function cθ( X::RealVec )
    tpfl=typeof(X[1])
    return X[4]/norm(X[2:4])
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   SINE OF PHI
#-----------------------------------------------------------------------
function sθ( X::RealVec )
    return X[3]/norm(X[3:4])
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   THE GORDON METRIC
#-----------------------------------------------------------------------
function gGordon( X::RealVec , n::Function=nIR 
                  , gfunc::Function=giso , U::RealVec=Float64[1;0;0;0] )

    tpfl=typeof(X[1])

    gs = gfunc( X )

    Unsq = abs(transpose(U)*gs*U)

    UU = uugen( U , X )

    return gs + ( 1-1/(n(X)^2) )*UU/Unsq
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#       END     Gordon.jl
#-----------------------------------------------------------------------