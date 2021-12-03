#-----------------------------------------------------------------------
module seval     # EVALUATION MODULE
#-----------------------------------------------------------------------

using LinearAlgebra , Statistics, Combinatorics, ForwardDiff
using DoubleFloats

include("type.jl")
include("metrics/WeakFieldIso.jl")
include("metrics/Minkowski.jl")
include("srl/FHC21.jl")
include("srl/RTC21.jl")
include("srl/mloc.jl")
include("outlier.jl")
include("geosol.jl")

Δψ0 = (π/180)*10    # Default angle margin

#-----------------------------------------------------------------------
#       STRUCTS
#-----------------------------------------------------------------------

"""
The TestCases datatype may be populated by the associated function of
the form:

    TestCases( par::NTuple , N::Int , np::Int , X::Array{RealMtx,1} , Xtar::Array{RealVec,1} )

The variables are defined in the following way:

    par     ::  NTuple              # Parameter tuple
    N       ::  Int                 # Number of generated test cases
    np      ::  Int                 # Number of emission points
    X       ::  Array{RealMtx,1}    # Emission points
    Xtar    ::  Array{RealVec,1}    # Target point

"""
struct TestCases                                    # TestCases datatype
    par     ::  NTuple              # Parameter tuple
    N       ::  Int                 # Number of cases
    np      ::  Int                 # Number of emission points
    X       ::  Array{RealMtx,1}    # Emission points
    Xtar    ::  Array{RealVec,1}    # Target point
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
The TestData datatype may be populated by the associated function of the
form:

    TestData(   par::NTuple , N::Int    , 
                X::Array{RealMtx,1}     ,   Xtar::Array{RealVec,1}     ,
                Xc::Array{RealVec,1}    ,   Xsc::Array{RealVec,1}      ,
                erh::RealVec    ,   erv::RealVec    ,   err::RealVec   ,
                erhC::RealVec   ,   ervC::RealVec   ,   errC::RealVec  ,
                Xc2::Array{RealVec,1}   ,   Xsc2::Array{RealVec,1}     ,
                erh2::RealVec   ,   erv2::RealVec   ,   err2::RealVec  ,
                erhC2::RealVec  ,   ervC2::RealVec  ,   errC2::RealVec )

The quantities are defined in the following way:

    par     ::  NTuple              # Parameter tuple
    N       ::  Int                 # Number of cases
    X       ::  Array{RealMtx,1}    # Emission points
    Xtar    ::  Array{RealVec,1}    # Target point
    Xc      ::  Array{RealVec,1}    # Flat spacetime position
    Xsc     ::  Array{RealVec,1}    # Squirrel position
    erh     ::  RealVec             # Horizontal error relative to Xsc
    erv     ::  RealVec             # Vertical error relative to Xsc
    err     ::  RealVec             # Total error relative to Xsc
    erhC    ::  RealVec             # Horizontal error relative to Xc
    ervC    ::  RealVec             # Vertical error relative to Xc
    errC    ::  RealVec             # Total error relative to Xc
    Xc2     ::  Array{RealVec,1}    # Flat spacetime position     (Aux.)
    Xsc2    ::  Array{RealVec,1}    # Squirrel position      (Auxiliary)
    erh2    ::  RealVec             # Horizontal errors      (Auxiliary)
    erv2    ::  RealVec             # Vertical error         (Auxiliary)
    err2    ::  RealVec             # Total error            (Auxiliary)
    erhC2   ::  RealVec             # Horizontal errors      (Auxiliary)
    ervC2   ::  RealVec             # Vertical error         (Auxiliary)
    errC2   ::  RealVec             # Total error            (Auxiliary)

The quantities suffixed with `2` are auxiliary quantities which are
needed in the four-point case, in which the special relativistic
location algorithms generally suffer from the bifurcation problem.

"""
struct TestData                                      # TestData datatype
    par     ::  NTuple              # Parameter tuple
    N       ::  Int                 # Number of cases
    X       ::  Array{RealMtx,1}    # Emission points
    Xtar    ::  Array{RealVec,1}    # Target point
    Xc      ::  Array{RealVec,1}    # Flat spacetime position
    Xsc     ::  Array{RealVec,1}    # Squirrel position
    erh     ::  RealVec             # Horizontal errors
    erv     ::  RealVec             # Vertical error
    err     ::  RealVec             # Total error
    erhC    ::  RealVec             # Horizontal errors
    ervC    ::  RealVec             # Vertical error
    errC    ::  RealVec             # Total error
    Xc2     ::  Array{RealVec,1}    # Flat spacetime position     (Aux.)
    Xsc2    ::  Array{RealVec,1}    # Squirrel position      (Auxiliary)
    erh2    ::  RealVec             # Horizontal errors      (Auxiliary)
    erv2    ::  RealVec             # Vertical error         (Auxiliary)
    err2    ::  RealVec             # Total error            (Auxiliary)
    erhC2   ::  RealVec             # Horizontal errors      (Auxiliary)
    ervC2   ::  RealVec             # Vertical error         (Auxiliary)
    errC2   ::  RealVec             # Total error            (Auxiliary)
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    tcfl( tc::TestCases , tpfl::DataType=Float64 )

The `tcfl` function changes the floating point type of `tc` to `tpfl`.

"""
function tcfl( tc::TestCases , tpfl::DataType=Float64 )
    par2    = tpfl.(tc.par)
    N       = tc.N
    np      = tc.np
    X2      = [zeros(tpfl,4,np) for _ in 1:N]
    Xtar2   = [zeros(tpfl,4) for _ in 1:N]

    for i=1:N 
        X2[i]     = tpfl.(tc.X[i])
        Xtar2[i]  = tpfl.(tc.Xtar[i])
    end

    return TestCases( par2 , N , np , X2 , Xtar2 )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#
#       GENERAL FUNCTIONS
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    roteuler( tpfl::DataType=Float64 )

The `roteuler` function constructs a random rotation matrix `R` and its 
inverse as a tuple `(R,inv(R))`.

"""
function roteuler( tpfl::DataType=Float64 )
    α = 2*pi*rand(tpfl)     # Euler angles
    β = 2*pi*rand(tpfl)
    γ = 2*pi*rand(tpfl)

    ( c1 , s1 ) = ( cos(α) , sin(α) )
    ( c2 , s2 ) = ( cos(β) , sin(β) )
    ( c3 , s3 ) = ( cos(γ) , sin(γ) )

    R = [     c2        -c3*s2           s2*s3          ;
            c1*s2   c1*c2*c3-s1*s3  -c3*s1-c1*c2*s3     ;
            s1*s2   c1*s3+c2*c3*s1   c1*c3-c2*s1*s3     ]

    return (R,inv(R))
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    angvec( u::RealVec , v::RealVec )

The `angvec` function computes the angle between the vectors `u` and 
`v`.

"""
function angvec( u::RealVec , v::RealVec )
    return acos( dot(u,v)/(norm(u)*norm(v)) )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    angmar( θ::Real , Δψ::Real )

The `angmar` function checks whether the angle `θ` falls in the range 
`0<θ<π/2 - |Δψ|`.

"""
function angmar( θ::Real , Δψ::Real )
    if θ <= π/2 && θ >+ 0
        ϕ = abs( π/2 - θ )
        if ϕ > Δψ
            return true
        else
            return false
        end
    else
        return false
    end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    slchk( X::RealMtx )

The `slchk` function checks whether the emission points in the matrix `X`
are spacelike separated with respect to the Minkowski metric `η_{μν}`.

"""
function slchk( X::RealMtx )
    tpfl = typeof(X[1])
    np = size(X)[2]

    if np == 1
        return true
    elseif np > 1
        Z = [zeros(tpfl,4) for _ in 1:np]

        for i=1:np
            Z[i] = X[:,i]
        end

        W = collect(combinations(Z,2))      # Combinations
        nc = length(W)
        Y = zeros(tpfl,nc)

        b = true
        for i=1:nc
            Y[i] = ηdot( W[i][1] - W[i][2] , W[i][1] - W[i][2] )
            if Y[i] > zero(tpfl)
                b = b * true
            else
                b = b * false
            end
        end

        return b
    end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#
#       POINT GENERATION
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#       Generates single random 3-component directional vector
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
"""
    vrgen( scale::Real , Δψ::Real , u::RealVec=Float64[1;0;0] )

The `vrgen` function generates a vector of length `scale` which makes
an angle `θ<π/2-|ΔΨ|`.

"""
function vrgen( scale::Real , Δψ::Real , u::RealVec=Float64[1;0;0] )
    # scale determines length, Δψ determines angle margins
    tpfl=typeof(scale)
    
    b = false
    
    v = zeros(tpfl,3)

    while b == false
        vx = 2*(rand(tpfl)-1/2)
        vy = 2*(rand(tpfl)-1/2)
        vz = 2*(rand(tpfl)-1/2)
        vl = sqrt( vx^2 + vy^2 + vz^2 )
        if vl>1
            b = false
        else
            b = true
        end
        v = tpfl[ vx ; vy ; vz ]/vl
        b = angmar(angvec(u,v),Δψ)
    end
    return scale*v
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    λiRscalc( xtar::RealVec , V::RealVec , Re::Real , Rs::Real )

The `λiRscalc` function computes the value of the affine parameter `λ`
for the test cases. Given the straight line `Y=Re*xtar+λ*V` passing 
through the point `Re*xtar`, it computes the value of `λ` at which the
spatial components of the line `Y[2:4]` have a norm `norm(Y[2:4])=Rs`.

"""
function λiRscalc( xtar::RealVec , V::RealVec , Re::Real , Rs::Real )
return  (
          -(Re*(V[1]*xtar[1] + V[2]*xtar[2] + V[3]*xtar[3])) +
          sqrt(abs(
            4*Rs^2*(V[1]^2 + V[2]^2 + V[3]^2) +
            Re^2*(
                  4*(V[1]*xtar[1] + V[2]*xtar[2] + V[3]*xtar[3])^2 -
                  4*(V[1]^2 + V[2]^2 + V[3]^2)*(xtar[1]^2 +
                  xtar[2]^2 + xtar[3]^2)
                 )
            ))/2
        )/(V[1]^2 + V[2]^2 + V[3]^2)
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    x3gen(tpfl::DataType=Float64)

The `x3gen` function generates a random 3-component position vector of
length `rell(x)` is the `WGS84` Earth reference ellipsoid radius.

"""
function x3gen(tpfl::DataType=Float64)
    x   = tpfl[ 1 ; 0 ; 0 ]
    rot = roteuler(tpfl)
    x = rot[1]*x
    return x*rell(x)
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    tidc( x::RealVec , v::RealVec , λ::Real , gfunc::Function )

The `tidc` function generates initial data for null geodesics given a
position vector `x`, velocity vector `v`, affine parameter `λ`, and a
metric function `gfunc`.

"""
function tidc( x::RealVec , v::RealVec , λ::Real , gfunc::Function )
    tpfl=typeof(x[1])

    X = zeros(tpfl,4)
    V = zeros(tpfl,4)

    X[2:4] = x
    V[2:4] = v

    V[1] = -λ

    V = nullenforcerp( V , X , gfunc )

    return vcat(X,V)
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#
#       PROBLEM GENERATOR
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    pgen( Rs::Real , gfunc::Function , tol::Real , np::Int=4 , 
          Δψ::Real=Δψ0 )

The `pgen` function generates a target point `Xtar` and `np` emission
points in a `4×np` matrix `X` such that the points in `X` lie on the
past light cone of `Xtar` with respect to the metric `gfunc`, and have
spatial radii values of `~Rs` (defined with respect to spatial origin).
The parameter `tol` is the tolerance parameter for the integration. This
function returns the tuple `(X,Xtar)`.

"""
function pgen( Rs::Real , gfunc::Function , tol::Real , np::Int=4 , 
               Δψ::Real=Δψ0 )
    tpfl=typeof(Rs)

    Z   = zeros(tpfl,8,np)
    Zi  = zeros(tpfl,8,np)

    X   = zeros(tpfl,4,np)

    Xc  = zeros(tpfl,4,4)

    xi  = x3gen(tpfl)
    vi  = zeros(tpfl,3,np)
    λi  = zeros(tpfl,np)

    for i=1:np
        B = true
        while B
            vi[:,i] = vrgen( one(tpfl) , Δψ , xi )
            λi[i]   = λiRscalc(xi/rell(xi),vi[:,i],rell(xi),Rs)
            vi[:,i] = λi[i]*vi[:,i]

            Zi[:,i] = tidc( xi , vi[:,i] , λi[i] , gfunc )
            Z[:,i]  = solveZ( Zi[:,i] , gfunc , tol , tol 
                              , AutoVern9(Rodas5()) )
            X[:,i]  = Z[1:4,i]

            # Check spacelike separation
            if i==1
                B = false
            elseif i>1 && i<4
                if slchk( X[:,1:i] )
                    B = false
                end
            elseif i>=4
                if slchk(X[:,1:i])
                    B = false
                end
            end
        end
    end

    return ( X , Zi[1:4,1] )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#
#       TEST CASE GENERATOR
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    gen( N::Int , g::Function , np::Int=4 , Δψ::Real=Δψ0 
              , tpfl::DataType=Float64 
              , REs ::  Real=1.4365276211950395e9       # Earth radius
              , RR  ::  Real=1.4365277e9    # Just above Earth surface
              , Rs  ::  Real=6e9            # Satellite radius
              , tolh::  Real=1e-14          # Tolerance for ODE solvers
              , tol ::  Real=1e-10          # Tolerance for ODE solvers
            )

The `gen` function generates `N` sets of target points `Xtar` and 
emission points `X` with respect to the metric `gfunc`. This function
returns a quantity of the TestCases datatype.

"""
function gen( N::Int , g::Function , np::Int=4 , Δψ::Real=Δψ0 
              , tpfl::DataType=Float64 
              , REs ::  Real=1.4365276211950395e9       # Earth radius
              , RR  ::  Real=1.4365277e9    # Just above Earth surface
              , Rs  ::  Real=6e9            # Satellite radius
              , tolh::  Real=1e-14          # Tolerance for ODE solvers
              , tol ::  Real=1e-10          # Tolerance for ODE solvers
            )

par = (REs,RR,Rs,tolh,tol)

tc = TestCases(  par , N , np , [zeros(tpfl,4,np) for _ in 1:N]    
                        , [zeros(tpfl,4) for _ in 1:N] )

for i=1:N 
    ZZ          = pgen( Rs , g , tolh , np , Δψ )
    tc.X[i]     = ZZ[1]
    tc.Xtar[i]  = ZZ[2]
    print("\r$i")
end

return tc
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#
#       EVALUATION
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    erH( Xtar::RealVec , ΔX::RealVec )

The `erH` function projects the vector `ΔX` in the direction orthogonal
to the vector `Xtar`. This is used to compute the horizontal error in
the terrestrial positioning problem.

"""
function erH( Xtar::RealVec , ΔX::RealVec )
    tpfl = typeof(Xtar[1])
    x    = copy(Xtar[2:4])
    Δx   = copy(ΔX[2:4])
    rv   = x/norm(x)
    δ    = one(tpfl)*(I(3))
    P    = δ - rv*transpose(rv)
    ΔE   = P*Δx

    return  norm(ΔE)
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    erV( Xtar::RealVec , ΔX::RealVec )

The `erV` function projects the vector `ΔX` in the direction along the
vector `Xtar`. This is used to compute the vertical error in the
terrestrial positioning problem.

"""
function erV( Xtar::RealVec , ΔX::RealVec )
    tpfl = typeof(Xtar[1])
    x    = copy(Xtar[2:4])
    Δx   = copy(ΔX[2:4])
    rv   = x/norm(x)

    return  dot(rv,Δx)
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    main(  tc::TestCases , sloc::Function , g::Function , Nx::Int=-1 , 
           tpfl::DataType=Double64 , tol::Real=1e-10 , ξ::Real=2e1 , 
           nb::Int=24 , ne::Real=6 )

This is the main evaluation function. It takes as input a collection of
test cases contained in a quantity of datatype `TestCases`, and
evaluates the squirrel locator function `sloc` on each test case
relative to a given metric function `g`. The function `main` returns
results in a quantity of type `TestData`.

The variable `Nx` is a limiter for the number of test cases to run. The
variable `tpfl` specifies the floating point precision for the
calculations. variables `tol`, `ξ`, `nb`,  and `ne` are inputs for the
`sloc` function, and correspond respectively to the integration
tolerance, outlier detection threshold, number of steps for the Broyden
algorithm, and number of emission points.

"""
function main(  tc::TestCases , sloc::Function , g::Function , 
                Nx::Int=-1 , tpfl::DataType=Float64 , tol::Real=1e-10 , 
                ξ::Real=2e1 , nb::Int=24 , ne::Real=6 )
tc2     = tcfl(tc,tpfl)
par     = tc2.par
REs     = par[1]
RR      = par[2]
Rs      = par[3]
tolh    = par[4]

if Nx==-1
N   = tc2.N
else
N   = Nx
end
np      = tc2.np

nc      = binomial(np,4)

X4      = zeros(tpfl,4,4)
ΔXsq    = zeros(tpfl,4)
ΔXcr    = zeros(tpfl,4)

ΔXs2    = [[zeros(tpfl,4)];[zeros(tpfl,4)]]
ΔXc2    = [[zeros(tpfl,4)];[zeros(tpfl,4)]]

d = size(tc2.X[1])[2]

td = TestData( par , N , tc2.X , tc2.Xtar 
         , [zeros(tpfl,4) for _ in 1:N]
         , [zeros(tpfl,4) for _ in 1:N]
         , zeros(tpfl,N) , zeros(tpfl,N) , zeros(tpfl,N)
         , zeros(tpfl,N) , zeros(tpfl,N) , zeros(tpfl,N)
         , [zeros(tpfl,4) for _ in 1:N]
         , [zeros(tpfl,4) for _ in 1:N]
         , zeros(tpfl,N) , zeros(tpfl,N) , zeros(tpfl,N)
         , zeros(tpfl,N) , zeros(tpfl,N) , zeros(tpfl,N)
         )

if d >= 5 && ne >= 5
    for i=1:N
        td.Xc[i]    = mlocator( td.X[i] )
        td.Xsc[i]   = sloc( td.X[i] , g , tol , nb , Double64 , ξ , ne )
        ΔXsq = td.Xtar[i] - td.Xsc[i]
        ΔXcr = td.Xtar[i] - td.Xc[i]

        td.erh[i]   = abs( erH( td.Xtar[i] , ΔXsq ) )
        td.erv[i]   = abs( erV( td.Xtar[i] , ΔXsq ) )
        td.err[i]   = norm( ΔXsq )

        td.erhC[i]  = abs( erH( td.Xtar[i] , ΔXcr ) )
        td.ervC[i]  = abs( erV( td.Xtar[i] , ΔXcr ) )
        td.errC[i]  = norm( ΔXcr )
        print("\r$i")
    end
elseif d == 4 || (d >= 4 && ne == 4)
    for i=1:N
        (td.Xc[i] ,td.Xc2[i])   = locator4FHC21( td.X[i] )
        (td.Xsc[i],td.Xsc2[i])  = sloc( td.X[i] , g , tol , nb , 
                                        Double64 , ξ , ne )
        ΔXsq  = td.Xtar[i] - td.Xsc[i]
        ΔXcr  = td.Xtar[i] - td.Xc[i]
        ΔXsq2 = td.Xtar[i] - td.Xsc2[i]
        ΔXcr2 = td.Xtar[i] - td.Xc2[i]

        td.erh[i]   = abs( erH( td.Xtar[i] , ΔXsq ) )
        td.erv[i]   = abs( erV( td.Xtar[i] , ΔXsq ) )
        td.err[i]   = norm( ΔXsq )

        td.erh2[i]   = abs( erH( td.Xtar[i] , ΔXsq2 ) )
        td.erv2[i]   = abs( erV( td.Xtar[i] , ΔXsq2 ) )
        td.er2r[i]   = norm( ΔXsq2 )

        td.erhC[i]  = abs( erH( td.Xtar[i] , ΔXcr ) )
        td.ervC[i]  = abs( erV( td.Xtar[i] , ΔXcr ) )
        td.errC[i]  = norm( ΔXcr )

        td.erhC2[i]  = abs( erH( td.Xtar[i] , ΔXcr2 ) )
        td.ervC2[i]  = abs( erV( td.Xtar[i] , ΔXcr2 ) )
        td.errC2[i]  = norm( ΔXcr2 )
        print("\r$i")
    end
end

return td
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
end     # END EVALUATION MODULE
#-----------------------------------------------------------------------