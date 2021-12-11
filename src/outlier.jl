#-----------------------------------------------------------------------
#
#   OUTLIER DETECTION FUNCTIONS
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    combX( X::RealMtx )

The `combX` function constructs and returns a vector of ``4×4`` matrices
consisting of all choices of 4 column vectors from the ``4×n_{\\rm e}`` matrix
`X`.

"""
function combX( X::RealMtx )
    tpfl=typeof(X[1,1])
    np = size(X)[2]

    if np==4
        return [X]
    elseif np>4
        Z = [zeros(tpfl,4) for _ in 1:np]   # Create containers

        for i=1:np
            Z[i] = copy(X[:,i])
        end

        W = collect(combinations(Z,4))      # Combinations
        nc = length(W)
        Y = [zeros(tpfl,4,4) for _ in 1:nc]

        for i=1:nc
            Y[i][:,1] = W[i][1]
            Y[i][:,2] = W[i][2]
            Y[i][:,3] = W[i][3]
            Y[i][:,4] = W[i][4]
        end

        return Y
    else
        return [zeros(tpfl,4,4)]
    end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    odetc( XC::Vector , ξ::Real=1e1 )

The `odetc` function takes a vector of vectors as input and returns a
vector of the vectors that differ by less than a threshold of `ξ`.

"""
function odetc( XC::Vector , ξ::Real=1e1 )
    tpfl    = typeof(XC[1][1])
    ns      = length(XC)

    Xc  = zeros(tpfl,4,ns)
    ΔX  = zeros(tpfl,4,ns)
    Xa  = zeros(tpfl,4)
    V   = zeros(tpfl,ns)
    BV  = zeros(Bool,ns)
    σ   = zeros(tpfl,ns)

    for i=1:ns
        Xc[:,i] = XC[i][:]
    end

    for i=1:4
        Xa[i] = median(Xc[i,:])  
    end

    for i=1:ns
        for j=1:4
            ΔX[j,i] = (Xc[j,i] - Xa[j])
        end
    end

    for i=1:ns
        V[i] = dot( ΔX[:,i] , ΔX[:,i] )
    end

    IV1 = partialsortperm( V , 1:ns )
    nr = 0

    if ns > 1
        for i=1:ns
            σ[i]    = sqrt(sum(V[IV1][1:i])/i)
            if i==1
                BV[i]   =  true
                nr = 1
            elseif i>1
                if σ[i] < ξ    
                    BV[i]   =  true
                    nr += 1
                end
            end
        end
    end
    if nr > 1
        IV2 = partialsortperm( BV , 1:ns , rev=true )
        return (XC[IV1[IV2]][1:nr],IV1[IV2])
    else
        return (XC[IV1][1:1],nr*ones(Int,ns))
    end
end     #---------------------------------------------------------------
