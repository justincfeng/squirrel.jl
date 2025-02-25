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

export ηdot, mnorm, ημν, δμν

# Default Earth parameters for isotropic metric
# [GM, J2, a]
const EARTH_ISO_PARAMS = [1.0, 1.0826300e-3, 1.438127773656399e9]

# Default Earth parameters for Kerr-Schild metric
# [a, GM]
const EARTH_KS_PARAMS = [738.0, 1.0]

# Default metric functions with Earth parameters
g = (X,p=EARTH_ISO_PARAMS,tpfl=Float64) -> gGordon(X,p)  # Gordon metric with Earth parameters
gks_earth = (X,p=EARTH_KS_PARAMS,tpfl=Float64) -> ge(X,p)        # Kerr-Schild metric with Earth parameters
giso_earth = (X,p=EARTH_ISO_PARAMS,tpfl=Float64) -> ge_iso(X,p)   # Isotropic metric with Earth parameters

# Parameterized metric functions
δ1=0.001
δ2=0.01
Patm=h->1.0
Pion=h->1.0

p0 = vcat(EARTH_ISO_PARAMS, [δ1, δ2])

gp  = (X, p=p0 , Pat=Patm , Pio=Pion ) -> gGordon(X, p[1:3], x->nIR(x,(h,θ,ϕ)->Δntp(h,θ,ϕ,p[4],p[5],Pat,Pio)), g)
gpc = (X, p=p0) -> gGordon(X, p[1:3], x->nIR(x,Δntpc))

#-----------------------------------------------------------------------
end # END MODULE METRIC
#-----------------------------------------------------------------------
