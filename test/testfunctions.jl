#-----------------------------------------------------------------------
#   MINKOWSKI METRIC
#-----------------------------------------------------------------------
η   = x->ημν()

#-----------------------------------------------------------------------
#   INITIAL DATA GENERATION
#-----------------------------------------------------------------------

function idg( ne::Int , Xc::RealVec=rand(4)/100 )
    V = rand(4,ne)
    V[1,:] = -ones(ne)

    for i=1:ne
        V[2:4,i] = V[2:4,i]/sqrt(dot(V[2:4,i],V[2:4,i]))
    end
    for i=1:ne
        @test ηdot(V[:,i],V[:,i]) ≈ 0 atol=1e-15
    end

    Y   = zeros(4,ne)
    Zi   = zeros(8,ne)

    for i=1:ne
        Y[:,i] = V[:,i] + Xc
    end
    for i=1:ne
        Zi[:,i] = vcat( Y[:,i] , Xc - Y[:,i] )
    end

    return (Y,Xc,Zi,V)
end
