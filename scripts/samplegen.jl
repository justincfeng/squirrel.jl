#-----------------------------------------------------------------------
#	SAMPLE GENERATION
#-----------------------------------------------------------------------

using LinearAlgebra, Serialization, BenchmarkTools

include("../src/squirrel.jl")
include("../src/metric.jl")

g  	= metric.g
gk 	= metric.ge

Nsamp   = 1000

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

tctloc	= dir*pfx*"-"*Nfs*sufx

Serialization.serialize(tctloc,squirrel.seval.tc2tup(tc))

#-----------------------------------------------------------------------
#	SAMPLES IN KERR GEOMETRY
#-----------------------------------------------------------------------

tck 	= squirrel.seval.gen(Nsamp,gk,6)

pfx	= "tck"

tckloc	= dir*pfx*"-"*Nfs*sufx

Serialization.serialize(tckloc,squirrel.seval.tc2tup(tck))

#-----------------------------------------------------------------------
