#-----------------------------------------------------------------------
#       BEGIN   atm.jl
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#   STANDARD ATMOSPHERIC PROFILE
#-----------------------------------------------------------------------
"""
    ΔnatmStd( h::Real )

The function `ΔnatmStd` is a profile for the atmospheric index of 
of refraction obtained from the U. S. standard atmosphere data. The
argument `h` is height from the Earth's surface in km.

"""
function ΔnatmStd( h::Real ) # h is height from surface in km
    tpfl=typeof(h)
    return  -222.66559677864387/
            (
                99.06212001632495 + 
                0.15782316597789148*(7.1540968904797175 + h)^2
            ) + 
            253.49880024119224/
            (
                112.75665221238798 + 
                0.17966910753975862*(7.156537366807337 + h)^2
            )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   REVISED EDLEN FORMULA
#-----------------------------------------------------------------------
"""
    ΔnREdlen( P::Real , T::Real , λ::Real )

The function `ΔnREdlen` provides the (revised) Edlen formula for 
computing the index of refraction from atmospheric parameters `P`, `T`,
and `λ`, corresponding to the respective pressure (in Pascals), 
temperature (in Celsius), and optical wavelength (in μm).

"""
function ΔnREdlen( P::Real , T::Real , λ::Real )
    tpfl=typeof(P)

    # dispersion relation
    disp = (tpfl(1e-8))*
           tpfl(8342.54 + 2406147.0/(130-λ^(-2)) + 15998.0/(38.9-λ^(-2)))

    return tpfl(P*disp/96095.43)*
           tpfl(
            ( 1. + (1e-8)*(0.601 - 0.00972*T)*P )
            /
            ( 1. + 0.0036610*T )
           )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#       END     atm.jl
#-----------------------------------------------------------------------