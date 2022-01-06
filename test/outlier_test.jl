#-----------------------------------------------------------------------
#   COMBINATION TEST
#-----------------------------------------------------------------------

for i=4:8
    @test length(combX( rand(4,i) ))    ≈   binomial(i,4)
    @test size(combX(rand(4,i))[1])     == size(rand(4,4))
    @test typeof(combX(rand(4,i))[1])   == typeof(rand(4,4))
end

#-----------------------------------------------------------------------
#   OUTLIER DETECTION TEST
#-----------------------------------------------------------------------

V = [2*(rand(4)/2-ones(4)) for _ in 1:10]
V[2] = (10^8)*2*(rand(4)/2-ones(4))
V[7] = (10^4)*ones(4)/2
V[8] = (10^6)*2*(rand(4)/2-ones(4))

OTL = odetc(V)

Vs  = V[OTL[2]]

ξ=1e1

@test OTL[2][8:10] == [7;8;2]

for i=1:7
    @test norm(Vs[i]) < norm(Vs[8])/ξ
    @test norm(Vs[i]) < norm(Vs[9])/ξ
    @test norm(Vs[i]) < norm(Vs[10])/ξ
end
