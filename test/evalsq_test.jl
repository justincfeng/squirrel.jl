#   INCOMPATIBLE WITH CURRENT VERSION

#-----------------------------------------------------------------------
#       roteuler & angvec TEST
#-----------------------------------------------------------------------

Rot = evalsq.roteuler()

v3 = rand(3)

@test (Rot[2]*(Rot[1]*v3)) ≈ v3 atol=1e-14

@test evalsq.angvec( Rot[1]*[1.;0.;0.]  , Rot[1]*[0.;1.;0.]  ) ≈ π/2
@test evalsq.angvec( Rot[1]*[1.;0.;0.]  , Rot[1]*[0.;0.;-1.] ) ≈ π/2
@test evalsq.angvec( Rot[1]*[0.;-1.;0.] , Rot[1]*[0.;0.;1.]  ) ≈ π/2
@test evalsq.angvec( Rot[1]*[1.;0.;0.]  , Rot[1]*[1.;1.;0.]  ) ≈ π/4

#-----------------------------------------------------------------------
#       angmar TEST
#-----------------------------------------------------------------------
Δψ0 = (π/180)*10

@test evalsq.angmar( 0 , Δψ0 )         == false
@test evalsq.angmar( π/2 , Δψ0 )       == false
@test evalsq.angmar( π/2-Δψ0/2 , Δψ0 ) == false

@test evalsq.angmar( Δψ0/2 , Δψ0 )
@test evalsq.angmar( Δψ0 , Δψ0 )
@test evalsq.angmar( 2*Δψ0 , Δψ0 )

@test evalsq.angmar( π/2-Δψ0 , Δψ0 )
@test evalsq.angmar( π/2-2*Δψ0 , Δψ0 )

#-----------------------------------------------------------------------
#       slchk TEST
#-----------------------------------------------------------------------

YE = idg(4)[1]

@test evalsq.slchk(YE)

YE[:,4] = YE[:,4] + [-10^4;0.;0.;0.]

@test evalsq.slchk(YE) == false

#-----------------------------------------------------------------------
#       TVchks TEST
#-----------------------------------------------------------------------

YE[:,1] = [0.;1.;0.;1.]
YE[:,2] = [0.;1.;0.;-1.]
YE[:,3] = [0.;-1.;0.;0.]
YE[:,4] = [0.;0.;1.;0.]

@test evalsq.TVchks( YE , 1e-9 )

YE[:,1] = [0.;1.;0.;1.]
YE[:,2] = [0.;1.;0.;-1.]
YE[:,3] = [0.;-1.;0.;0.]
YE[:,4] = [0.;0.;0.;0.]

@test evalsq.TVchks( YE , 1e-9 ) == false

#-----------------------------------------------------------------------
#       cereal checker TEST
#-----------------------------------------------------------------------

YE = idg(4,zeros(4))[1]

@test evalsq.crlchks( YE )

YE[:,4] = YE[:,4] + [-10^4;0.;0.;0.]

@test evalsq.crlchks( YE )  == false

#-----------------------------------------------------------------------
#       vrgen TEST
#-----------------------------------------------------------------------

u = rand(3)
sc = 100.
v = evalsq.vrgen( sc , Δψ0 , u )

@test sqrt(dot(v,v)) ≈ sc atol=5e-14
@test evalsq.angmar(evalsq.angvec(u,v),Δψ0)

#-----------------------------------------------------------------------
#       irscalc TEST
#-----------------------------------------------------------------------

v = u/sqrt(dot(u,u))

u = rand(3)
x = u/sqrt(dot(u,u))

@test evalsq.irscalc( x , v , 0  , 100 ) ≈ 100 atol=5e-14
@test evalsq.irscalc( x , x , 10 , 100 ) ≈ 90  atol=5e-14
@test evalsq.irscalc( x , x , 50 , 100 ) ≈ 50  atol=5e-14

#-----------------------------------------------------------------------
#       tidc TEST
#-----------------------------------------------------------------------

λ = evalsq.irscalc( x , v , 10 , 100 )
v = λ*v
Zi = evalsq.tidc( x , v , λ , η )

@test Zi[2:4] ≈ x
@test Zi[6:8] ≈ v

#-----------------------------------------------------------------------
#       pgen TEST
#-----------------------------------------------------------------------
ne    = 6
ntest = 5

P = evalsq.pgen( 6e9 , η , 1e-14 , ne , Δψ0 , true , 3 , ntest ) ;

@test (P[1]-P[3][2])./P[1] ≈ zeros(4,ne) atol=1e-13
for i=3:ntest
    @test (P[3][i-1]-P[3][i])./P[3][i-1] ≈ zeros(4,ne) atol=1e-13
end
@test (P[1]-P[3][ntest])./P[1] ≈ zeros(4,ne) atol=1e-13

#-----------------------------------------------------------------------
#       gen & main TEST
#-----------------------------------------------------------------------

nev = 2

nb	= 24
tol	= 1e-9
ξ1	= 1e4
ξ2	= 1e1

for k=1:10

tc = evalsq.gen(3,η,ne)
typeof(tc) == evalsq.TestCases

td = evalsq.main(tc,squirrel.mlocator,η,nb,tol,nev,Float64,true,ξ1,ξ2,ne)

@test length(td.erh) ≈ nev

Xer = [zeros(4) for _ in 1:nev]

for i=1:nev
    Xer[i][2:4] = (td.Xsc[i][2:4]-td.Xtar[i][2:4])./td.Xtar[i][2:4]
    @test Xer[i] ≈ zeros(4) atol=tol
end

end
