using Combinatorics

const RealVec{T<:Real} = Array{T,1}     ;
const RealMtx{T<:Real} = Array{T,2}     ;

function breakapart( X::RealMtx )  # Break apart function
    tpfl = typeof(X[1,1])
    d  = size(X)[1]
    np = size(X)[2]
    Z = [zeros(tpfl,d) for _ in 1:np]

    for i=1:np
        Z[i] = X[:,i]
    end

    return Z
end

X  = rand(4,4)

X4 = breakapart(X)

w = collect( combinations( X4 , 3 ) )

w = w[[4,3,2,1]]    # This reverses the order of elements
