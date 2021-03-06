#-----------------------------------------------------------------------
#   TOLERANCES
#-----------------------------------------------------------------------

δtol = 1e-8
tol  = 1e-13

#-----------------------------------------------------------------------
#   GEODESIC SOLVER TEST
#-----------------------------------------------------------------------

Z0 = [0.;0.;0.;0.;1.;1.;0.;0.]      # Index raised in p's here
ZF = [1.;1.;0.;0.;-1.;1.;0.;0.]     # Index lowered in p's here
@test squirrel.gsolve(  Z0[1:4] , Z0[5:8] , η , δtol 
                        , AutoVern7(Rodas5()) ) ≈ ZF

Z0 = [0.;0.;0.;0.;1.;0.;1.;0.]      # Index raised in p's here
ZF = [1.;0.;1.;0.;-1.;0.;1.;0.]     # Index lowered in p's here
@test squirrel.gsolve(  Z0[1:4] , Z0[5:8] , η , δtol 
                        , AutoVern7(Rodas5()) ) ≈ ZF

Z0 = [0.;0.;0.;0.;1.;0.;0.;1.]      # Index raised in p's here
ZF = [1.;0.;0.;1.;-1.;0.;0.;1.]     # Index lowered in p's here
@test squirrel.gsolve(  Z0[1:4] , Z0[5:8] , η , δtol 
                        , AutoVern7(Rodas5()) ) ≈ ZF

#-----------------------------------------------------------------------
#   V34 TEST
#-----------------------------------------------------------------------

v3 = rand(3)
@test length(squirrel.V34(v3)) ≈ 4
@test squirrel.V34(v3)[2:4] == v3

#-----------------------------------------------------------------------
#   VidF TEST
#-----------------------------------------------------------------------

Zi = zeros(8,4)
Zi[6:8,1] = ones(3)
Zi[6:8,2] = 2*ones(3)
Zi[6:8,3] = 3*ones(3)
Zi[6:8,4] = 4*ones(3)

Vid = squirrel.VidF( Zi )

@test Vid[1:3] ≈ ones(3)
@test Vid[4:6] ≈ 2*ones(3)
@test Vid[7:9] ≈ 3*ones(3)
@test Vid[10:12] ≈ 4*ones(3)

#-----------------------------------------------------------------------
#   INITIAL DATA GENERATION
#-----------------------------------------------------------------------

X   = [-0.9938120402728282 -0.9938120402728282 -0.9938120402728282 -0.9938120402728282; 
        0.43329524598249924 0.6468681369385507 0.615089937176196 0.5122363212193108; 
        0.738110796898415 0.7597257376287113 0.2839892724192345 0.0952119165738411; 
        0.5369302954880543 0.15569297625663295 0.7506786499110514 0.8654592211384085]

Xc  = [0.006187959727171788, 0.008549046454540871, 0.004581697263911051, 0.006340316336591607]

Zi  = [-0.9938120402728282 -0.9938120402728282 -0.9938120402728282 -0.9938120402728282; 
        0.43329524598249924 0.6468681369385507 0.615089937176196 0.5122363212193108; 
        0.738110796898415 0.7597257376287113 0.2839892724192345 0.0952119165738411; 
        0.5369302954880543 0.15569297625663295 0.7506786499110514 0.8654592211384085; 
        1.0 1.0 1.0 1.0; 
       -0.4247461995279584 -0.6383190904840098 -0.6065408907216551 -0.5036872747647699; 
       -0.733529099634504 -0.7551440403648003 -0.27940757515532344 -0.09063021930993004; 
       -0.5305899791514627 -0.14935265992004135 -0.7443383335744598 -0.8591189048018169]

V   = [-1.0 -1.0 -1.0 -1.0; 
        0.4247461995279584 0.6383190904840098 0.6065408907216551 0.5036872747647699; 
        0.733529099634504 0.7551440403648003 0.27940757515532344 0.09063021930993004; 
        0.5305899791514627 0.14935265992004135 0.7443383335744598 0.8591189048018169]

Vid  = vcat(-V[2:4,1],-V[2:4,2],-V[2:4,3],-V[2:4,4])

#-----------------------------------------------------------------------
#   ZERO FUNCTION TEST
#-----------------------------------------------------------------------

@test squirrel.zF( Vid , Zi , η , δtol ) ≈ zeros(12) atol=tol

#-----------------------------------------------------------------------
#   GEJAC TEST
#-----------------------------------------------------------------------

for i=1:4
    local res = squirrel.gejac( Zi[1:4,i] , Zi[5:8,i] , η , δtol )
    @test res[1][1:4] ≈ Xc atol=tol
    @test res[1][6:8] ≈ -V[2:4,i] atol=tol
    @test res[2][1,:] ≈ -V[2:4,i] atol=tol
end

#-----------------------------------------------------------------------
#   JACOBIAN CALCULATOR TEST
#-----------------------------------------------------------------------
global resx = squirrel.geocJ( Zi , η , δtol )

@test resx[1][6:8,:]       ≈ -V[2:4,:] atol=tol

for i=1:4
    @test resx[1][1:4,i]   ≈ Xc atol=tol
end

@test resx[2][1,1:3]       ≈ -V[2:4,1] atol=tol
@test resx[2][5,1:3]       ≈ -V[2:4,1] atol=tol
@test resx[2][9,1:3]       ≈ -V[2:4,1] atol=tol

@test resx[2][1,4:6]       ≈ V[2:4,2] atol=tol
@test resx[2][5,7:9]       ≈ V[2:4,3] atol=tol
@test resx[2][9,10:12]     ≈ V[2:4,4] atol=tol

@test resx[2][5:12,  4:6]  ≈ zeros(8,3)
@test resx[2][1:4,   7:9]  ≈ zeros(4,3)
@test resx[2][9:12,  7:9]  ≈ zeros(4,3)
@test resx[2][1:8, 10:12]  ≈ zeros(8,3)

@test resx[2][2:4,1:3]     ≈ 1.0*I(3)
@test resx[2][6:8,1:3]     ≈ 1.0*I(3)
@test resx[2][10:12,1:3]   ≈ 1.0*I(3)
@test resx[2][2:4,4:6]     ≈ -1.0*I(3)
@test resx[2][6:8,7:9]     ≈ -1.0*I(3)

#-----------------------------------------------------------------------
#   zFc TEST
#-----------------------------------------------------------------------

Zf = ones(8,4)
Zf[:,2] = -ones(8)
Zf[:,3] = -2*ones(8)
Zf[:,4] = -3*ones(8)

@test squirrel.zFc( Zf )[1:4]  ≈ 2*ones(4)
@test squirrel.zFc( Zf )[5:8]  ≈ 3*ones(4)
@test squirrel.zFc( Zf )[9:12] ≈ 4*ones(4)

#-----------------------------------------------------------------------
#   INITIAL DATA CORRECTOR TEST
#-----------------------------------------------------------------------

Zf = squirrel.idf( Zi , η , δtol , 24 ) 
@test Zf[1:4] ≈ Zi[1:4] atol=δtol*mean(Zi)/10
@test Zf[6:8] ≈ Zi[6:8] atol=δtol*mean(Zi)/10

#-----------------------------------------------------------------------
#   INITIAL DATA CONSTRUCTOR
#-----------------------------------------------------------------------

gk   = x->gks(x,0.5) ;

XE   = [4.8552418244084965 4.845669922167221 4.856689375629473 4.8206419799059494 4.823208150136851 4.856202508967677  ; 
        7.299483910135271 6.8054911167005425 7.1810230696314425 6.153377751048404 7.746460952899268 7.227281564130913  ; 
        10.349124252528748 9.710191586980727 10.434190218681847 10.14678773174398 9.282597162251408 10.217530411064208 ; 
        5.265994805337563 6.471968677435302 5.6150675825992185 5.951777443574898 5.391088065951989 6.0625298198886775  ] ;
      
Xtar = [5.855953201888279 ; 6.78444938070708 ; 9.525600927112148 ; 5.492794659774428] ;

#-----------------------------------------------------------------------
#   SINGLE LOCATOR TEST
#-----------------------------------------------------------------------

Xs = squirrel.locator4( XE[:,1:4] , 
                        squirrel.locator4FHC22(XE[:,1:4])[1] , gk ,
                        1e-10 , 24 , false )

@test Xs ≈ Xtar  atol=5e-14

#-----------------------------------------------------------------------
#   MULTI LOCATOR TEST
#-----------------------------------------------------------------------

ne  = 6

Xsc  = squirrel.locator(  XE , gk , 1e-10 , ne , 24 , 2e1 , Double64  )

@test Xsc ≈ Xtar  atol=5e-14
