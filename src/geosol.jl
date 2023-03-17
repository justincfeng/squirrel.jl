#-----------------------------------------------------------------------
#       BEGIN   geosol.jl
#-----------------------------------------------------------------------

using OrdinaryDiffEq

#-----------------------------------------------------------------------
#       INITIAL DATA NULL ENFORCER FUNCTIONS
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    nullenforcerf( V0::RealVec , X::RealVec , gfunc::Function ) 

The `nullenforcerf` function takes a vector `V0` and modifies the time
or first component so that it is null at the point `X` with respect to
`gfunc` and future directed; the function returns the modified vector.

"""
function nullenforcerf( V0::RealVec , X::RealVec , gfunc::Function )  
    tpfl=typeof(V0[1])
    nx = length(V0)

    Vn = copy(V0)

    g = gfunc(X)

    if g[1,1] < zero(tpfl)
        Vn[1] =  -((g[1,2]*Vn[2] + g[1,3]*Vn[3] + g[1,4]*Vn[4] + 
        sqrt( 
              4*(g[1,2]*Vn[2] + g[1,3]*Vn[3] + g[1,4]*Vn[4])^2 - 
              4*g[1,1]*(g[2,2]*Vn[2]^2 + 2*g[2,3]*Vn[2]*Vn[3] + 
              g[3,3]*Vn[3]^2 + 2*g[2,4]*Vn[2]*Vn[4] + 
              2*g[3,4]*Vn[3]*Vn[4] + g[4,4]*Vn[4]^2)
            )/2)/g[1,1])

        return Vn

    elseif g[1,1] == zero(tpfl)
        Vn[1] =  (-(g[2,2]*Vn[2]^2) - 2*g[2,3]*Vn[2]*Vn[3] - 
                    g[3,3]*Vn[3]^2 - 2*g[2,4]*Vn[2]*Vn[4] - 
                    2*g[3,4]*Vn[3]*Vn[4] - 
                    g[4,4]*Vn[4]^2)/(2*(g[1,2]*Vn[2] + 
                    g[1,3]*Vn[3] + g[1,4]*Vn[4]))

        return Vn
    else 
        print("Need g_{11} to be negative or zero.")
        return zeros(tpfl,2*nx)
    end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    nullenforcerp( V0::RealVec , X::RealVec , gfunc::Function ) 

The `nullenforcerp` function takes a vector `V0` and modifies the time
or first component so that it is null at the point `X` with respect to
`gfunc` and past directed; the function returns the modified vector.

"""
function nullenforcerp( V0::RealVec , X::RealVec , gfunc::Function )
    tpfl=typeof(V0[1])
    nx = length(V0)

    Vn = copy(V0)

    g = gfunc(X)

    if g[1,1] < zero(tpfl)
        Vn[1] =  (-2*g[1,2]*Vn[2] - 2*g[1,3]*Vn[3] - 2*g[1,4]*Vn[4] + 
        sqrt(4*(g[1,2]*Vn[2] + g[1,3]*Vn[3] + g[1,4]*Vn[4])^2 - 
        4*g[1,1]*(g[2,2]*Vn[2]^2 + 2*g[2,3]*Vn[2]*Vn[3] + 
        g[3,3]*Vn[3]^2 + 2*g[2,4]*Vn[2]*Vn[4] + 2*g[3,4]*Vn[3]*Vn[4] + 
        g[4,4]*Vn[4]^2)))/(2*g[1,1])

        return Vn

    elseif g[1,1] == zero(tpfl)
        Vn[1] =  (-(g[2,2]*Vn[2]^2) - 2*g[2,3]*Vn[2]*Vn[3] - 
                    g[3,3]*Vn[3]^2 - 2*g[2,4]*Vn[2]*Vn[4] - 
                    2*g[3,4]*Vn[3]*Vn[4] - 
                    g[4,4]*Vn[4]^2)/(2*(g[1,2]*Vn[2] + 
                    g[1,3]*Vn[3] + g[1,4]*Vn[4]))

        return Vn
    else 
        print("Need g_{11} to be negative or zero.")
        return zeros(tpfl,2*nx)
    end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    HamGeo( Z::RealVec , gfunc::Function )

The `HamGeo` function computes the geodesic Hamiltonian 
``H = \\frac{1}{2} g^{μν} p_μ p_ν``
from the metric ``g_{μν}`` provided in the function `gfunc` and the 
spacetime position `x` and momentum `p` encoded in `Z=(x,p)`.

"""
function HamGeo( Z::RealVec , gfunc::Function )
    tpfl=typeof(Z[1])

    ( x , p ) = ( Z[1:4] , Z[5:8] )

    gin = inv(gfunc(x))

    return tpfl(1/2)*(transpose(p)*gin*p)
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    Jsympl( Zarg::RealVec )

The `Jsympl` function effectively applies the symplectic matrix
``J^{αβ}`` to the vector `Zarg` (the vector ``{∂H}/{∂z^β}``).

"""
function Jsympl( Zarg::RealVec )
    tpfl=typeof(Zarg[1])
    n2 = length(Zarg)
    Z = zeros(tpfl,n2)

    if iseven(n2)
        n = Int(n2 / 2)
        for i=1:n
            Z[i]    = Zarg[n+i] 
            Z[n+i]  = - Zarg[i] 
        end

        return Z
    else
        return Z
    end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    ZdotGeo( Z::RealVec , gfunc::Function )

The `ZdotGeo` function computes the gradient of the Hamiltonian and
applies the symplectic operator to the gradient.

"""
function ZdotGeo( Z::RealVec , gfunc::Function )
    tpfl=typeof(Z[1])

    VH = ForwardDiff.gradient(z->HamGeo(z,gfunc),Z)           
    return Jsympl( VH )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    ZF!( dZ , Z , gfunc::Function )

The function `ZF` converts the function `ZdotGeo` into a function with
mutable arguments.

"""
function ZF!( dZ , Z , gfunc::Function )
    copyto!(dZ, ZdotGeo( Z, gfunc ))
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    solveZ( Z0::RealVec , gfunc::Function , δ1::Real , δ2::Real ,
            integrator=AutoVern7(Rodas5()) , δt=0 )

The function `solveZ` performs the integration of geodesics given
the initial data `Z0`, and the metric function `gfunc`. The integration
is performed using the integrator specified by the variable `integrator`
using the relative tolerance parameter `δ1` and the absolute tolerance 
parameter `δ2`. The parameter `δt` determines the timestep parameter 
`dt` in ODEProblem. 

The function `solveZ` outputs only the endpoint of the solution for the
geodesic equation.

"""
function solveZ( Z0::RealVec , gfunc::Function , δ1::Real , δ2::Real , 
    integrator=AutoVern7(Rodas5()) , δt=0 ; 
    cb=CallbackSet() )
tpfl=typeof(Z0[1])

tspan = (0.0,1.0) 

Zo = vcat(Z0[1:4],gfunc(Z0[1:4])*Z0[5:8]) # Lowering index of p 

F = (dZ,Z,p,t)->ZF!( dZ , Z , gfunc )

if δt == 0
slv = solve( ODEProblem(F,Zo,tspan),integrator,reltol=δ1
    ,abstol=δ2,save_everystep=false,callback=cb)
else
slv = solve( ODEProblem(F,Zo,tspan),integrator,reltol=δ1
    ,abstol=δ2,dt=δt,save_everystep=false,callback=cb )
end

return last(slv.u)
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#       END     geosol.jl
#-----------------------------------------------------------------------