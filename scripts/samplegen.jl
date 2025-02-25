#-----------------------------------------------------------------------
#	SAMPLE GENERATION
#-----------------------------------------------------------------------

using LinearAlgebra, Serialization, BenchmarkTools

include("../src/squirrel.jl")
include("../src/metric.jl")

# Define parameter vectors
p_iso = metric.EARTH_ISO_PARAMS
p_ks = metric.EARTH_KS_PARAMS

# Define metric functions with parameters
g = (x,p=p_iso)->metric.g(x,p)
gk = (x,p=p_ks)->metric.gks(x,p)

Nsamp   = 100

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

tc	= squirrel.seval.gen(Nsamp,g,p_iso,6)

pfx	= "tct"

tctloc	= dir*pfx*"-"*Nfs*sufx

Serialization.serialize(tctloc,squirrel.seval.tc2tup(tc))

#-----------------------------------------------------------------------
#	SAMPLES IN KERR GEOMETRY
#-----------------------------------------------------------------------

tck 	= squirrel.seval.gen(Nsamp,gk,p_ks,6)

pfx	= "tck"

tckloc	= dir*pfx*"-"*Nfs*sufx

Serialization.serialize(tckloc,squirrel.seval.tc2tup(tck))

#-----------------------------------------------------------------------
