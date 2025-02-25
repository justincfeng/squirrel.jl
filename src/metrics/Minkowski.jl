#-----------------------------------------------------------------------
#       BEGIN   Minkowski.jl
#-----------------------------------------------------------------------

export ηdot, mnorm, ημν, δμν

#-----------------------------------------------------------------------
#   THE MINKOWSKI PRODUCT (RECTANGULAR COORDINATES)
#-----------------------------------------------------------------------
"""
    ηdot( V1::RealVec , V2::RealVec )

The function `ηdot` takes two vectors ``V_1`` and ``V_2`` of equal
dimension, and computes their Minkowski product ``η(V_1,V_2)``, with the
assumption that the first element of the vectors corresponds to the time
component.

"""
function ηdot( V1::RealVec , V2::RealVec )  # Computes Minkowski product
    nv1 = length(V1)
    nv2 = length(V2)

    met = -V1[1]*V2[1]

    if nv1==nv2
        for i=2:nv1
            met += V1[i]*V2[i]
        end
        return met
    else
        print("Vectors are of a different dimension.")
        return 0*V[1]
    end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   THE MINKOWSKI NORM
#-----------------------------------------------------------------------
"""
    mnorm( V ) 

The function `mnorm` computes the Minkowski norm, which is
mathematically equivalent to the absolute value of the square root of
the Minkowski product ``√|η(V_1,V_2)|``. However, this function computes
the result according to the formula used in the hypot function.

"""
function mnorm( V )
    l = length(V)

    t = abs(V[1])
    x = norm(V[2:l])

    if t > x
        return abs(t - ( x/t )*( x / ( 1 + √( 1 - (x/t)^2 ) ) ))
    elseif x > t
        return abs(x - ( t/x )*( t / ( 1 + √( 1 - (t/x)^2 ) ) ))
    else
        return zero( typeof(V[1]) )
    end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   THE MINKOWSKI METRIC (RECTANGULAR COORDINATES)
#-----------------------------------------------------------------------
"""
    ημν( tpfl::DataType=Float64 , p=Float64[] , dim::Int=4 )

The function `ημν` constructs the components of the Minkowski metric.
The first argument `tpfl` specifies the floating point datatype
(typically Float64 or Double64 if one uses the DoubleFloats package),
the second argument `dim` specifies the dimension, and the third argument
`p` is an optional parameter vector (included for API consistency but not used
since Minkowski metric is parameter-free). The default values
are `tpfl=Float64`, `dim=4`, and `p=nothing`:

    julia> ημν() == ημν( Float64 , 4 , nothing )
        true

"""
function ημν( tpfl::DataType=Float64 , p=Float64[] , dim::Int=4 )
    gη = one(tpfl)*(I(dim))
    gη[1,1] = -gη[1,1]
    return Matrix(gη)
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   KRONECKER DELTA 
#-----------------------------------------------------------------------
"""
    δμν( tpfl::DataType=Float64 , p=Float64[] , dim::Int=4 )

The function `δμν` returns the components of an identity matrix of 
dimension `dim` using the floating-point datatype `tpfl`. The parameter
vector `p` is included for API consistency but not used since the identity
matrix is parameter-free.
"""
function δμν( tpfl::DataType=Float64 , p=Float64[] , dim::Int=4)
    return one(tpfl)*(I(dim))
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#       END     Minkowski.jl
#-----------------------------------------------------------------------
