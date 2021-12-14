#-----------------------------------------------------------------------
#       BEGIN   atmios.jl
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#   PSEUDO-EPSTEIN FUNCTION
#-----------------------------------------------------------------------
"""
    pEp( h , hc , B , Nmax )

The function `pEp` implements the pseudo-Epstein line shape function
used to approximate the Epstein function used to construct the 
ionospheric electron density profiles.

"""
function pEp( h , hc , B , Nmax )
    return  (Nmax/16.0)*(   1/(1 + ((h - hc)/(2.0*B))^2) + 
                            1/(1 + ((h - hc)/(4.0*B))^2)^2 + 
                            1/(1 + ((h - hc)/(6.0*B))^2)^3 + 
                            1/(1 + ((h - hc)/(7.0*B))^2)^4
                        )^2
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   LINE SHAPE FUNCTION
#-----------------------------------------------------------------------
"""
    LSF( x , xo , σ )

The function `LSF` implements a line shape function employed in the
perturbation models; it has a faster falloff than the Lorentzian.

"""
function LSF( x , xo , σ )
    return  ( σ^2 / ( σ^2 + (x-xo)^2 ) )*
            ( σ^4 / ( σ^4 + (x-xo)^4 ) )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   INDEX OF REFRACTION FOR IONOSPHERE
#-----------------------------------------------------------------------
"""
    Δnios( h::Real , θ::Real , ϕ::Real )

The function `Δnios` provides a crude model for the ionospheric
electron density profile. The variable `h` represents height from the
Earth's surface in units of km.

"""
function Δnios( h::Real , θ::Real , ϕ::Real )
    # This is a crude model which emulates ionospheric effects
    # h is height from surface in km
    tpfl=typeof(h)
    return  ( 
            pEp(h,300,50,4.02411e-5) +    # F layer
            pEp(h,130,30,1.00603e-5) +    # E layer
            pEp(h,75,5,4.02411e-6)        # D layer
            )  
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   TOTAL INDEX OF REFRACTION
#-----------------------------------------------------------------------
"""
    Δntot( h::Real , θ::Real , ϕ::Real )

The function `Δntot` is the total index of refraction profile, which
is a sum of atmospheric and ionospheric index of refraction profiles. 
The variable `h` represents height from the Earth's surface in units of 
km.

"""
function Δntot( h::Real , θ::Real , ϕ::Real )
    # h is height from surface in km
    tpfl=typeof(h)
    return  ΔnatmStd( h ) + Δnios( h , θ , ϕ )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   PERTURBATION MODELS
#-----------------------------------------------------------------------
"""
    P( h::Real , h0::RealVec , σ::RealVec )

The function `P` constructs the perturbation model, which consists of
a sum of line shape functions `LSF`, with the parameters provided in
the vector arguments `h0` and `σ`.

"""
function P( h::Real , h0::RealVec , σ::RealVec )
    tpfl=typeof(h)
    n = length(h0)
    return sum(LSF(h,h0[i],σ[i]) for i=1:n)
end

#-----------------------------------------------------------------------
#   TOTAL INDEX OF REFRACTION WITH PERTURBATION
#-----------------------------------------------------------------------
"""
    Δntp( h::Real , θ::Real , ϕ::Real , δ1::Real=0.001
               , δ2::Real=0.010 , Patm::Function=h->1.0
               , Pion::Function=h->1.0 )

The function `Δntp` adds the perturbations provided in the perturbation
model to the atmospheric and ionospheric profiles defined earlier.

"""
function Δntp( h::Real , θ::Real , ϕ::Real , δ1::Real=0.001
               , δ2::Real=0.010 , Patm::Function=h->1.0
               , Pion::Function=h->1.0 )
    # h is height from surface in km
    tpfl=typeof(h)

    return  ΔnatmStd( h ) * ( 1 + δ1*Patm(h) ) + 
            Δnios(h,θ,ϕ) * ( 1 + δ2*Pion(h) )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   COMPUTE INDEX OF REFRACTION (SPHERICAL COORDINATES)
#-----------------------------------------------------------------------
"""
    nIRs( RE::Real , r::Real , θ::Real , ϕ::Real , Δnf::Function )

The function `nIRs` computes the index of refraction in spherical
coordinates from the index of refraction profile function `Δnf` in
spherical coordinates. 

"""
function nIRs( RE::Real , r::Real , θ::Real , ϕ::Real , Δnf::Function )
    # RE and r in units of Earth mass
    tpfl=typeof(RE)
    um = tpfl(4.435e-6)                     # Conversion factor to km
    h = (r - RE)*um                         # h in units of km
    return  one(tpfl) + Δnf(h,θ,ϕ)
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   COMPUTE INDEX OF REFRACTION
#-----------------------------------------------------------------------
"""
    nIR( X::RealVec , Δnf::Function=Δntot )

The function `nIR` computes the index of refraction in Cartesian
coordinates; the profile is anchored to the WGS84 ellipsoid defined by
the `rell` function.

"""
function nIR( X::RealVec , Δnf::Function=Δntot )
    # RE and r in units of Earth mass

    xs = SphericalFromCartesian()(X[2:4])

    return  nIRs( rell(X[2:4]) , norm(X[2:4]) , xs.θ , xs.ϕ 
                  , Δnf )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#       END     atmios.jl
#-----------------------------------------------------------------------