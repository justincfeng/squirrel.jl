#   INCOMPATIBLE WITH CURRENT VERSION

#-----------------------------------------------------------------------
#   MINKOWSKI PRODUCT TEST
#-----------------------------------------------------------------------

tv  = [1.;0.;0.;0.]
xv  = [0.;1.;0.;0.]
yv  = [0.;0.;1.;0.]
zv  = [0.;0.;0.;1.]

@test cereal.η(tv,tv) ≈ -1
@test cereal.η(xv,xv) ≈ 1
@test cereal.η(yv,yv) ≈ 1
@test cereal.η(zv,zv) ≈ 1

@test cereal.η(tv,xv) ≈ 0
@test cereal.η(tv,yv) ≈ 0
@test cereal.η(tv,zv) ≈ 0
@test cereal.η(xv,yv) ≈ 0
@test cereal.η(xv,zv) ≈ 0
@test cereal.η(yv,zv) ≈ 0

k1  = π^2*[1.;1.;0.;0.]
k2  = π^2*[1.;0.;1.;0.]
k3  = π^2*[1.;0.;0.;1.]

kn1 = π^2*[1.;-1.;0.;0.]
kn2 = π^2*[1.;0.;-1.;0.]
kn3 = π^2*[1.;0.;0.;-1.]

@test cereal.η(k1,k1) ≈ 0
@test cereal.η(k2,k2) ≈ 0
@test cereal.η(k3,k3) ≈ 0

@test cereal.η(kn1,kn1) ≈ 0
@test cereal.η(kn2,kn2) ≈ 0
@test cereal.η(kn3,kn3) ≈ 0

#-----------------------------------------------------------------------
#   FRAME TEST
#-----------------------------------------------------------------------

ne = 4
ID4 = idg(ne)
Y4 = ID4[1]
Y4[:,ne] = Y4[:,ne] + [-10^4;0.;0.;0.]

@test cereal.Frame(Y4,false) ≈ zeros(4,3)

#-----------------------------------------------------------------------
#   NORMVEC TEST
#-----------------------------------------------------------------------

ID4 = idg(4,zeros(4))
@test cereal.NormVecF(cereal.Frame(ID4[1])) ≈ [1,0,0,0]

#-----------------------------------------------------------------------
#   LORENTZ TRANSFORMATION TEST
#-----------------------------------------------------------------------

NV = [1,0,0,0] + rand(4)/20
NV = NV/sqrt(abs(ηdot(NV,NV)))
Λ = cereal.LTM(NV)

@test Λ*NV ≈ [1,0,0,0] atol=1e-13

#-----------------------------------------------------------------------
#   IPFINDER TEST
#-----------------------------------------------------------------------

@test cereal.IPfinder(ID4[1]) ≈ zeros(4) atol=1e-13

#-----------------------------------------------------------------------
#   LOCATOR TEST
#-----------------------------------------------------------------------

X = cerealtest.epgen()
@test cerealtest.single( 1e-13 , cereal.locator(X) , X )

for i=1:Int(100)
    X2 = cerealtest.epgen()
    @test cerealtest.single( 1e-13 , cereal.locator(X2) , X2 )
end
