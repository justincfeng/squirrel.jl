#-----------------------------------------------------------------------
#       DataType TEST
#-----------------------------------------------------------------------

tpfl    = Float64
par     = (1.4365276211950395e9,1.4365277e9,6e9,1e-14,1e-10)
N       = 10
np      = 5 
X       = [rand(tpfl,4,np) for _ in 1:N] 
Xtar    = [rand(tpfl,4) for _ in 1:N]

tcR =  squirrel.seval.TestCases(  par , N , np , X , Xtar )

@test tcR.par   == par
@test tcR.N     == N
@test tcR.np    == np
@test tcR.X     == X
@test tcR.Xtar  == Xtar

Xc      = [rand(tpfl,4) for _ in 1:N]
Xsc     = [rand(tpfl,4) for _ in 1:N]
erh     = rand(tpfl,N)
erv     = rand(tpfl,N)
err     = rand(tpfl,N)
erhC    = rand(tpfl,N)
ervC    = rand(tpfl,N)
errC    = rand(tpfl,N)
Xc2     = [rand(tpfl,4) for _ in 1:N]
Xsc2    = [rand(tpfl,4) for _ in 1:N]
erh2    = rand(tpfl,N)
erv2    = rand(tpfl,N)
err2    = rand(tpfl,N)
erhC2   = rand(tpfl,N)
ervC2   = rand(tpfl,N)
errC2   = rand(tpfl,N)

tdR =  squirrel.seval.TestData( par , N , X , Xtar , Xc , Xsc , erh , 
                erv , err , erhC , ervC , errC , Xc2 , Xsc2 , erh2 , 
                erv2 , err2 , erhC2 , ervC2 , errC2 )

@test tdR.par   == par
@test tdR.N     == N
@test tdR.X     == X
@test tdR.Xtar  == Xtar
@test tdR.Xc    == Xc
@test tdR.Xsc   == Xsc
@test tdR.erh   == erh
@test tdR.erv   == erv
@test tdR.err   == err
@test tdR.erhC  == erhC
@test tdR.ervC  == ervC
@test tdR.errC  == errC
@test tdR.Xc2   == Xc2
@test tdR.Xsc2  == Xsc2
@test tdR.erh2  == erh2
@test tdR.erv2  == erv2
@test tdR.err2  == err2
@test tdR.erhC2 == erhC2
@test tdR.ervC2 == ervC2
@test tdR.errC2 == errC2

#-----------------------------------------------------------------------
#       DataType converter TEST
#-----------------------------------------------------------------------

tcRD = squirrel.seval.tcfl( tcR , Double64 )

@test typeof(tcRD.X[1][1,1])  == Double64
@test typeof(tcRD.Xtar[1][1]) == Double64

tcT = squirrel.seval.tc2tup( tcR )

@test tcT == (par,N,np,X,Xtar)
@test tcT == (tcR.par,tcR.N,tcR.np,tcR.X,tcR.Xtar)

tdT = squirrel.seval.td2tup( tdR )

@test tdT == (par,N,X,Xtar,Xc,Xsc,erh,erv,err,erhC,ervC,errC,Xc2,Xsc2,
              erh2,erv2,err2,erhC2,ervC2,errC2)

@test tdT == (tdR.par,tdR.N,tdR.X,tdR.Xtar,tdR.Xc,tdR.Xsc,tdR.erh,
              tdR.erv,tdR.err,tdR.erhC,tdR.ervC,tdR.errC,tdR.Xc2,
              tdR.Xsc2,tdR.erh2,tdR.erv2,tdR.err2,tdR.erhC2,tdR.ervC2,
              tdR.errC2)

tcRb = squirrel.seval.tup2tc( tcT )

tdRb = squirrel.seval.tup2td( tdT )

@test tcRb == tcR

@test tdRb == tdR

#-----------------------------------------------------------------------
#       roteuler & angvec TEST
#-----------------------------------------------------------------------

Rot = squirrel.seval.roteuler()

v3 = rand(3)

@test (Rot[2]*(Rot[1]*v3)) ≈ v3 atol=1e-14

@test squirrel.seval.angvec( Rot[1]*[1.;0.;0.]  , 
                             Rot[1]*[0.;1.;0.]  ) ≈ π/2
@test squirrel.seval.angvec( Rot[1]*[1.;0.;0.]  , 
                             Rot[1]*[0.;0.;-1.] ) ≈ π/2
@test squirrel.seval.angvec( Rot[1]*[0.;-1.;0.] , 
                             Rot[1]*[0.;0.;1.]  ) ≈ π/2
@test squirrel.seval.angvec( Rot[1]*[1.;0.;0.]  , 
                             Rot[1]*[1.;1.;0.]  ) ≈ π/4

#-----------------------------------------------------------------------
#       angmar TEST
#-----------------------------------------------------------------------
Δψ0 = (π/180)*10

@test squirrel.seval.angmar( 0 , Δψ0 )         == false
@test squirrel.seval.angmar( π/2 , Δψ0 )       == false
@test squirrel.seval.angmar( π/2-Δψ0/2 , Δψ0 ) == false

@test squirrel.seval.angmar( Δψ0/2 , Δψ0 )
@test squirrel.seval.angmar( Δψ0 , Δψ0 )
@test squirrel.seval.angmar( 2*Δψ0 , Δψ0 )

@test squirrel.seval.angmar( π/2-Δψ0 , Δψ0 )
@test squirrel.seval.angmar( π/2-2*Δψ0 , Δψ0 )

#-----------------------------------------------------------------------
#       slchk TEST
#-----------------------------------------------------------------------

YE = idg(4)[1]

@test squirrel.seval.slchk(YE)

YE[:,4] = YE[:,4] + [-10^4;0.;0.;0.]

@test squirrel.seval.slchk(YE) == false

#-----------------------------------------------------------------------
#       vrgen TEST
#-----------------------------------------------------------------------

u = rand(3)
sc = 100.
v = squirrel.seval.vrgen( sc , Δψ0 , u )

@test sqrt(dot(v,v)) ≈ sc atol=5e-14
@test squirrel.seval.angmar(squirrel.seval.angvec(u,v),Δψ0)

#-----------------------------------------------------------------------
#       λiRscalc TEST
#-----------------------------------------------------------------------

v = u/sqrt(dot(u,u))

u = rand(3)
x = u/sqrt(dot(u,u))

@test squirrel.seval.λiRscalc( x , v , 0  , 100 ) ≈ 100 atol=5e-14
@test squirrel.seval.λiRscalc( x , x , 10 , 100 ) ≈ 90  atol=5e-14
@test squirrel.seval.λiRscalc( x , x , 50 , 100 ) ≈ 50  atol=5e-14

#-----------------------------------------------------------------------
#       tidc TEST
#-----------------------------------------------------------------------

λ = squirrel.seval.λiRscalc( x , v , 10 , 100 )
v = λ*v
Zi = squirrel.seval.tidc( x , v , λ , η )

@test Zi[2:4] ≈ x
@test Zi[6:8] ≈ v

#-----------------------------------------------------------------------
#       erv & erh TEST
#-----------------------------------------------------------------------

U =  [0.24981918733261255, 0.6256573311421523, 0.8453503802138926, 
      0.3093630156499423]

ΔU = [0.000981491116663281, 0.0008969460327706518, 
      0.00034983811644343675, 0.00011308545853025453]

@test   squirrel.seval.erV(U,ΔU) ≈ 
        dot(U[2:4],ΔU[2:4])/norm(U[2:4]) atol=tol

@test   √(squirrel.seval.erV(U,ΔU)^2 + squirrel.seval.erH(U,ΔU)^2) ≈ 
        norm(ΔU[2:4]) atol=tol

#-----------------------------------------------------------------------
#       pgen TEST
#-----------------------------------------------------------------------

ne    = 6
ntest = 5

P = squirrel.seval.pgen( 6e9 , η , 1e-14 , ne , Δψ0 )

@test size(P[1]) == (4,ne)
@test size(P[2]) == (4,)

#-----------------------------------------------------------------------
#       gen & main TEST
#-----------------------------------------------------------------------

nev = 2

nb	= 24
tol	= 1e-9
ξ	= 1e1

for k=1:10

tc = squirrel.seval.gen(3,η,ne)
typeof(tc) == squirrel.seval.TestCases

td = squirrel.seval.main(tc,squirrel.locator,η,nev,Float64,tol,ξ,nb,ne)

@test length(td.erh) ≈ nev
@test length(td.erv) ≈ nev
@test length(td.err) ≈ nev

Xer = [zeros(4) for _ in 1:nev]

for i=1:nev
    Xer[i][2:4] = (td.Xsc[i][2:4]-td.Xtar[i][2:4])./td.Xtar[i][2:4]
    @test Xer[i] ≈ zeros(4) atol=tol
end

end
