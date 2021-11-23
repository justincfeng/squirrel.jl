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

g   = X->gGordon(X)
gp  =   (X,δ1=0.001,δ2=0.01,Patm=h->1.0,Pion=h->1.0)->
        gGordon(X,x->nIR(x,(h,θ,ϕ)->Δntp(h,θ,ϕ,δ1,δ2)))
gpc = X->gGordon(X,x->nIR(x,Δntpc))

#-----------------------------------------------------------------------
end # END MODULE METRIC
#-----------------------------------------------------------------------
