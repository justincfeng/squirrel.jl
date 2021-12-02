#-----------------------------------------------------------------------
#	SAMPLE GENERATION
#-----------------------------------------------------------------------

using LinearAlgebra, Serialization, BenchmarkTools

include("../src/squirrel.jl")
include("../src/metric.jl")

g  	= metric.g
gk 	= metric.ge

Nsamp   = 100
# N	= 100000

nb	= 24
tol	= 1e-14
ξ1	= 1e-18
ξ2	= 1e1	
sufx 	= ""
dir 	= "../res/"
Nfs	= string(Nsamp)
	
#-----------------------------------------------------------------------
#	SAMPLES IN ANALOGUE GEOMETRY WITH ATMOSPHERIC & IONOSPHERIC EFFECTS
#-----------------------------------------------------------------------

tc	= squirrel.seval.gen(Nsamp,g,6)

pfx	= "tct"

tct 	= (Nsamp,tc.par,tc.X,tc.Xtar)

tctloc	= dir*pfx*"-"*Nfs*sufx

Serialization.serialize(tctloc,tct)

#-----------------------------------------------------------------------
#	SAMPLES IN KERR GEOMETRY
#-----------------------------------------------------------------------

tck 	= squirrel.seval.gen(Nsamp,gk,6)

pfx	= "tck"

tctk 	= (Nsamp,tck.par,tck.X,tck.Xtar)

tckloc	= dir*pfx*"-"*Nfs*sufx

Serialization.serialize(tckloc,tctk)

#-----------------------------------------------------------------------
