#-----------------------------------------------------------------------
#	SAMPLE GENERATION
#-----------------------------------------------------------------------

using LinearAlgebra, Serialization, BenchmarkTools
using CoordinateTransformations, Geodesy

include("../src/evalsq.jl")
include("../src/cereal.jl")
include("../src/squirrel.jl")
include("../src/geocoord.jl")
include("../src/metric.jl")

g  	= metric.g
gk 	= metric.ge

N	= 100
# N	= 100000

nb	= 24
tol	= 1e-14
ξ1	= 1e4
ξ2	= 1e1	
sufx 	= ""
dir 	= "../res/"
Ns	= string(N)
	
#-----------------------------------------------------------------------
#	SAMPLES IN ANALOGUE GEOMETRY WITH ATMOSPHERIC & IONOSPHERIC EFFECTS
#-----------------------------------------------------------------------

tc	= evalsq.gen(N,g,6)

pfx	= "tct"

tct 	= (N,tc.par,tc.X,tc.Xtar)

tctloc	= dir*pfx*"-"*Ns*sufx

Serialization.serialize(tctloc,tct)

#-----------------------------------------------------------------------
#	SAMPLES IN KERR GEOMETRY
#-----------------------------------------------------------------------

tck 	= evalsq.gen(N,gk,6)

pfx	= "tck"

tctk 	= (N,tck.par,tck.X,tck.Xtar)

tckloc	= dir*pfx*"-"*Ns*sufx

Serialization.serialize(tckloc,tctk)

#-----------------------------------------------------------------------
