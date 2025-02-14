#-----------------------------------------------------------------------
module metric
#-----------------------------------------------------------------------
using LinearAlgebra
using CoordinateTransformations
using LegendrePolynomials

include("type.jl")

include("atmios/atm.jl")
include("atmios/atmios.jl")

include("metrics/Minkowski.jl")
include("metrics/KerrSchild.jl")
include("metrics/WeakFieldIso.jl")
include("metrics/Gordon.jl")

# Default Earth parameters for isotropic metric
# [GM, J2, a]
const EARTH_ISO_PARAMS = [1.0, 1.0826300e-3, 1.438127773656399e9]

# Default Earth parameters for Kerr-Schild metric
# [a, GM]
const EARTH_KS_PARAMS = [738.0, 1.0]

# Default metric functions with Earth parameters
g   = X -> ge_gordon(X)  # Gordon metric with Earth parameters
gks = X -> ge(X)        # Kerr-Schild metric with Earth parameters
giso = X -> ge_iso(X)   # Isotropic metric with Earth parameters

# Parameterized metric functions
gp  = (X, δ1=0.001, δ2=0.01, Patm=h->1.0, Pion=h->1.0) ->
        gGordon(X, EARTH_ISO_PARAMS, x->nIR(x,(h,θ,ϕ)->Δntp(h,θ,ϕ,δ1,δ2)))
gpc = X -> gGordon(X, EARTH_ISO_PARAMS, x->nIR(x,Δntpc))

#-----------------------------------------------------------------------
end # END MODULE METRIC
#-----------------------------------------------------------------------
