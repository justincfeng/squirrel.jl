#-----------------------------------------------------------------------
"""
    multivec( X::RealMtx , k::Int )

This function takes a ``m×n`` matrix `X`, and for `k```<n``,
constructs a vector of ``m×```k` matrices constructed from all choices
of `k` columns from the ``m×n`` matrix `X`.

"""
function multivec( X::RealMtx , k::Int )  
    tpfl = typeof(X[1,1])
    d    = size(X)[1]
    np   = size(X)[2]
    Z    = [zeros(tpfl,d) for _ in 1:np]

    for i=1:np
        Z[i] = X[:,i]
    end

    w   = collect(combinations(Z,k))

    l   = length(w)
    z   = [zeros(tpfl,d,k) for _ in 1:l]
    for i=1:l, j=1:k
        z[i][:,j] = w[i][j]
    end

    return z
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    mlocator( X::RealMtx )

The function `mlocator` computes a single intersection point ``X_c`` 
from an ``m×n`` matrix `X` formed from `n>4` emission points, which form
the column vectors of `X`.

"""
function mlocator( X::RealMtx )
    #   Calculates location for more than five emission points
    tpfl = typeof(X[1,1])

    l = size(X)

    if l[2] == 5
        return locator5RTC21(X)
    elseif l[2] > 5
        W   = multivec(X,5)
        k   = length(W)

        w   = [zeros(tpfl,4) for _ =1:2*k ]

        for a=1:k
            w[a] = locator5RTC21(W[a])
        end

        #   Re-sort location points according to smallest Minkowski norm
        δW = zeros(tpfl,k)
        for a=1:k
            for i=1:4
                δW[a] += abs(ηdot(X[:,i]-w[a],X[:,i]-w[a]))
            end
        end
        iW = sortperm(δW)

        return w[iW][1]
    else 
        return zeros(tpfl,4)
    end
end     #--------------------------------------------------------------- 
