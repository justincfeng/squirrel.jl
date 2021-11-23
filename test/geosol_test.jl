#-----------------------------------------------------------------------
#   TEST NULL ENFORCER FUNCTIONS
#-----------------------------------------------------------------------

X   = rand(4)

v       = zeros(4)
v[2:4]  = rand(3)
vnorm  = ηdot(v,v)
nvf = nullenforcerf( v , X , η )
nvp = nullenforcerp( v , X , η )

@test nvf[1] > 0.
@test ηdot(nvf,nvf) ≈ 0 atol=1e-15

@test nvp[1] < 0.
@test ηdot(nvp,nvp) ≈ 0 atol=1e-15

#-----------------------------------------------------------------------
#   TEST HAMILTONIAN
#-----------------------------------------------------------------------

Z = [0.;0.;0.;0.;1.;0.;0.;0.]
@test HamGeo( Z , η ) ≈ -0.5

Z = [0.;0.;0.;0.;0.;1.;0.;0.]
@test HamGeo( Z , η ) ≈ 0.5

v3 = rand(3)
Z[6:8] = v3
@test HamGeo( Z , η ) ≈ dot(v3,v3)/2

#-----------------------------------------------------------------------
#   TEST SYMPLECTIC OPERATOR
#-----------------------------------------------------------------------

Z = ones(8)

@test Jsympl(zeros(8)) ≈ zeros(8)
@test Jsympl(π*Z)[1:4] ≈ π*ones(4)
@test Jsympl(π*Z)[5:8] ≈ -π*ones(4)

#-----------------------------------------------------------------------
#   TEST HAMILTON EQUATIONS
#-----------------------------------------------------------------------

@test ZdotGeo( π*Z , η )[1]     ≈ -π
@test ZdotGeo( π*Z , η )[2:4]   ≈ π*ones(3)
@test ZdotGeo( π*Z , η )[5:8]   ≈ zeros(4)

#-----------------------------------------------------------------------
#   TEST SOLVEZ FUNCTION
#-----------------------------------------------------------------------

Z0 = [0.;0.;0.;0.;1.;1.;0.;0.]      # Index raised in p's here
ZF = [1.;1.;0.;0.;-1.;1.;0.;0.]     # Index lowered in p's here
@test solveZ( Z0 , η , 1e-9 , 1e-9 , AutoVern7(Rodas5()) ) ≈ ZF

Z0 = [0.;0.;0.;0.;1.;0.;1.;0.]      # Index raised in p's here
ZF = [1.;0.;1.;0.;-1.;0.;1.;0.]     # Index lowered in p's here
@test solveZ( Z0 , η , 1e-9 , 1e-9 , AutoVern7(Rodas5()) ) ≈ ZF

Z0 = [0.;0.;0.;0.;1.;0.;0.;1.]      # Index raised in p's here
ZF = [1.;0.;0.;1.;-1.;0.;0.;1.]     # Index lowered in p's here
@test solveZ( Z0 , η , 1e-9 , 1e-9 , AutoVern7(Rodas5()) ) ≈ ZF
