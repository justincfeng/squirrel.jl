#-----------------------------------------------------------------------
#   CHECK MINKOWSKI METRIC FUNCTIONS USED HERE
#-----------------------------------------------------------------------

tv  = [1.;0.;0.;0.]
xv  = [0.;1.;0.;0.]
yv  = [0.;0.;1.;0.]
zv  = [0.;0.;0.;1.]

@test ηdot(tv,tv) ≈ -1
@test ηdot(xv,xv) ≈ 1
@test ηdot(yv,yv) ≈ 1
@test ηdot(zv,zv) ≈ 1

@test ηdot(tv,xv) ≈ 0
@test ηdot(tv,yv) ≈ 0
@test ηdot(tv,zv) ≈ 0
@test ηdot(xv,yv) ≈ 0
@test ηdot(xv,zv) ≈ 0
@test ηdot(yv,zv) ≈ 0

k1  = π^2*[1.;1.;0.;0.]
k2  = π^2*[1.;0.;1.;0.]
k3  = π^2*[1.;0.;0.;1.]

kn1 = π^2*[1.;-1.;0.;0.]
kn2 = π^2*[1.;0.;-1.;0.]
kn3 = π^2*[1.;0.;0.;-1.]

@test ηdot(k1,k1) ≈ 0
@test ηdot(k2,k2) ≈ 0
@test ηdot(k3,k3) ≈ 0

@test ηdot(kn1,kn1) ≈ 0
@test ηdot(kn2,kn2) ≈ 0
@test ηdot(kn3,kn3) ≈ 0

@test ημν()[1,1] ≈ -1
@test ημν()[2,2] ≈ 1
@test ημν()[3,3] ≈ 1
@test ημν()[4,4] ≈ 1

@test ημν()[1,2] ≈ 0
@test ημν()[1,3] ≈ 0
@test ημν()[1,4] ≈ 0
@test ημν()[2,3] ≈ 0
@test ημν()[2,4] ≈ 0
@test ημν()[3,4] ≈ 0
