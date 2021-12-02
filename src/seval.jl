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
# TestCases datatype

This datatype may be populated by the associated function of the form:

    TestCases( par::NTuple , N::Int , np::Int , X::Array{RealMtx,1} , Xtar::Array{RealVec,1} )

The variable 'par' is a tuple of real parameters, 'N' is the number of generated test cases, 'np' is the number of emission points, 'X' is the matrix of emission points, and 'Xtar' is the target point.

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
# TestData datatype

This datatype may be populated by the associated function of the form:

    TestData( par::NTuple , N::Int , X::Array{RealMtx,1} , Xtar::Array{RealVec,1} , Xc::Array{RealVec,1} , Xsc::Array{RealVec,1} , erh::RealVec , erv::RealVec , err::RealVec , erhC::RealVec , ervC::RealVec , errC::RealVec , Xc2::Array{RealVec,1} , Xsc2::Array{RealVec,1} , erh2::RealVec , erv2::RealVec , err2::RealVec , erhC2::RealVec , ervC2::RealVec , errC2::RealVec )

The variable 'par' is a tuple of real parameters, 'N' is the number of generated test cases, 'X' is the matrix of emission points, 'Xtar' is the target point, 'Xc' is the solution computed in flat spacetime, and 'Xsc' is the solution computed in curved spacetime. The quantities 'erh', 'erv'. and 'err' are the respective horizontal, vertical and total errors between 'Xsc' and 'Xtar', and the corresponding errors suffixed with 'C' are are constructed from 'Xc' and 'Xtar'. The quantities suffixed with '2' are auxiliary quantities which are needed in the four-point case.

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

#-----------------------------------------------------------------------
#
#       GENERAL FUNCTIONS
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#       Random rotation generator
#-----------------------------------------------------------------------
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
#       Angle between vectors
#-----------------------------------------------------------------------
function angvec( u::RealVec , v::RealVec )
    return acos( dot(u,v)/(norm(u)*norm(v)) )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#       Angle margin checker
#-----------------------------------------------------------------------
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
#
#       FUNCTIONS FOR CHECKING srl5 COMPATIBILITY
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#       Check spacelike separation
#-----------------------------------------------------------------------
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
#       Intersection with Rs calculator
#-----------------------------------------------------------------------
function irscalc( x::RealVec , dx::RealVec , Re::Real , Rs::Real )
    return  (
                -(Re*(dx[1]*x[1] + dx[2]*x[2] + dx[3]*x[3])) +
                sqrt(abs(
                    4*Rs^2*(dx[1]^2 + dx[2]^2 + dx[3]^2) +
                    Re^2*(
                           4*(dx[1]*x[1] + dx[2]*x[2] + dx[3]*x[3])^2 - 
                           4*(dx[1]^2 + dx[2]^2 + dx[3]^2)*(x[1]^2 + 
                           x[2]^2 + x[3]^2)
                         )
                    ))/2
            )/(dx[1]^2 + dx[2]^2 + dx[3]^2)
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#       Generates random 3-component position vector
#-----------------------------------------------------------------------
function x3gen(tpfl::DataType=Float64)
    x   = tpfl[ 1 ; 0 ; 0 ]
    rot = roteuler(tpfl)
    x = rot[1]*x
    return x*rell(x)
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#       Initial data construction
#-----------------------------------------------------------------------
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

function pgen( Rs::Real , gfunc::Function , δ::Real , n::Int=4 
               , Δψ::Real=Δψ0 , ctest::Bool=false , δf::Real=2 
               , ntest::Int=4 )
    tpfl=typeof(Rs)

    Z   = zeros(tpfl,8,n)
    Zi  = zeros(tpfl,8,n)

    X   = zeros(tpfl,4,n)

    Xc  = zeros(tpfl,4,4)

    xi  = x3gen(tpfl)
    vi  = zeros(tpfl,3,n)
    λi  = zeros(tpfl,n)

    for i=1:n
        B = true
        while B
            vi[:,i] = vrgen( one(tpfl) , Δψ , xi )
            λi[i]   = irscalc(xi/rell(xi),vi[:,i],rell(xi),Rs)
            vi[:,i] = λi[i]*vi[:,i]

            Zi[:,i] = tidc( xi , vi[:,i] , λi[i] , gfunc )
            Z[:,i]  = solveZ( Zi[:,i] , gfunc , δ , δ 
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

    if ctest
        Zn   = [zeros(tpfl,8,n) for _ in 1:ntest]
        Xn   = [zeros(tpfl,4,n) for _ in 1:ntest]

        for j=1:ntest
            δ = δ/δf
            for i=1:n
                Zn[j][:,i]  = solveZ( Zi[:,i] , gfunc , δ , δ 
                              , AutoVern9(Rodas5()) )
                Xn[j][:,i]  = Zn[j][1:4,i]
            end
        end
        return ( X , Zi[1:4,1] , Xn )
    else
        return ( X , Zi[1:4,1] )
    end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#
#       TEST CASE GENERATOR
#
#-----------------------------------------------------------------------

function gen( N::Int , g::Function , np::Int=4 , Δψ::Real=Δψ0 
              , tpfl::DataType=Float64 
              , REs ::  Real=1.4365276211950395e9       # Earth radius
              , RR  ::  Real=1.4365277e9    # Just above Earth surface
              , Rs  ::  Real=6e9            # Satellite radius
              , tolh::  Real=1e-14          # Tolerance for ODE solvers
              , tol ::  Real=1e-9           # Tolerance for ODE solvers
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
#       Horizontal Error
#-----------------------------------------------------------------------
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
#       Vertical Error
#-----------------------------------------------------------------------
function erV( Xtar::RealVec , ΔX::RealVec )
    tpfl = typeof(Xtar[1])
    x    = copy(Xtar[2:4])
    Δx   = copy(ΔX[2:4])
    rv   = x/norm(x)

    return  dot(rv,Δx)
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#       MAIN EVALUATION FUNCTION
#-----------------------------------------------------------------------
function main(  tc::TestCases , sloc::Function , g::Function ,
                nb::Int=24 , tol::Real=1e-10 , Nx::Int=-1 , 
                tpfl::DataType=Float64 , ξ::Real=2e1 , ne::Real=6 )
par     = tc.par
REs     = par[1]
RR      = par[2]
Rs      = par[3]
tolh    = par[4]

if Nx==-1
N   = tc.N
else
N   = Nx
end
np      = tc.np

nc      = binomial(np,4)

X4      = zeros(tpfl,4,4)
ΔXsq    = zeros(tpfl,4)
ΔXcr    = zeros(tpfl,4)

ΔXs2    = [[zeros(tpfl,4)];[zeros(tpfl,4)]]
ΔXc2    = [[zeros(tpfl,4)];[zeros(tpfl,4)]]

d = size(tc.X[1])[2]

td = TestData( par , N , tc.X , tc.Xtar 
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