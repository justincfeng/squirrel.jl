#-----------------------------------------------------------------------
#   COMBINATION TEST
#-----------------------------------------------------------------------

for i=4:8
    @test length(combX( rand(4,i) ))    ≈   binomial(i,4)
    @test size(combX(rand(4,i))[1])     == size(rand(4,4))
    @test typeof(combX(rand(4,i))[1])   == typeof(rand(4,4))
end

#-----------------------------------------------------------------------
#   CBT TEST
#-----------------------------------------------------------------------

Xc = [(rand(4),false) for _ in 1:20]
Xc[14:20] = [(rand(4),true) for _ in 1:7]

@test CBT(Xc,2) ≈ 7

#-----------------------------------------------------------------------
#   OUTLIER DETECTION TEST
#-----------------------------------------------------------------------

V = [2*(rand(4)/2-ones(4)) for _ in 1:10]
V[2] = (10^8)*2*(rand(4)/2-ones(4))
V[7] = (10^4)*ones(4)/2
V[8] = (10^6)*2*(rand(4)/2-ones(4))

OTL = odetc(V)

@test OTL[2][8:10] == [7;8;2]

for i=1:7
    @test dot(OTL[1][i],OTL[1][i])      <  dot(OTL[1][8],OTL[1][8])
end

#-----------------------------------------------------------------------
#   CEREAL ERROR FILTER TEST
#-----------------------------------------------------------------------

ne = 4
ID4 = idg(ne)
Y4 = ID4[1]
Y4[:,ne] = Y4[:,ne] + [-10^4;0.;0.;0.]

w4 = cefilt( Y4 ) ;
@test w4[1][1] ≈ zeros(4,4) atol=1e-15

ne = 6
ID6 = idg(ne)
Y6 = ID6[1]
Y6[:,ne] = Y6[:,ne] + [-10^4;0.;0.;0.]

w6 = cefilt( Y6 ) ;
@test length(w6[1]) ≈ 5

ne = 7
ID7 = idg(ne)
Y7 = ID7[1]
Y7[:,ne] = Y7[:,ne] + [-10^4;0.;0.;0.]

w7 = cefilt( Y7 ) ;
@test length(w7[1]) ≈ 15

#-----------------------------------------------------------------------
#   CEREAL OUTLIER FILTER
#-----------------------------------------------------------------------

@test length(cofilt(Y4)[1]) ≈ 1
@test length(cofilt(Y6)[1]) ≈ 5
@test length(cofilt(Y7)[1]) ≈ 15
