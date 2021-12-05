#-----------------------------------------------------------------------
#		EVALUATION SCRIPT
#           It is recommended that you run this script with 4 threads
#           enabled.
#-----------------------------------------------------------------------

using LinearAlgebra, Serialization, BenchmarkTools

include("../src/squirrel.jl")
include("../src/metric.jl")

g  	= metric.g          # Gordon metric with standard parameters
gk 	= metric.ge         # Kerr-Schild with Earth parameters

Neval	= 1000            # Number of test cases to evaluate
Nsamp	= 100000        # Number of test cases in generated sample file

nb      = 24            # Number of steps for Broyden solver
tol	    = 1e-10         # Tolerance for ODE solver (OrdinaryDiffEq.jl)
tpfl    = Float64       # Floating point precision to use
ξ	    = 2e1           # Threshold for outlier detection code

sufx 	= ""                # Filename string suffix
dir 	= "../res/"         # Results directory
Nes	    = string(Neval)     # Turn Neval into a string
Nfs	    = string(Nsamp)     # Turn Nfile into a string
	
#-----------------------------------------------------------------------
#		DESERIALIZE GENERATED SAMPLE FILES (ANALOGUE GEOMETRY)
#-----------------------------------------------------------------------
	
pfx	= "tct" ;
tctloc	= dir*pfx*"-"*Nfs*sufx ;

tct	= Serialization.deserialize(tctloc) ;
	
tc	= squirrel.seval.tup2tc( tct ) ;

#-----------------------------------------------------------------------
#		DESERIALIZE GENERATED SAMPLE FILES (KERR GEOMETRY)
#-----------------------------------------------------------------------
	
pfx	= "tck" ;
tckloc	= dir*pfx*"-"*Nfs*sufx ;

tctk	= Serialization.deserialize(tckloc) ;
	
tck	 = squirrel.seval.tup2tc( tctk ) ;

#-----------------------------------------------------------------------
ne	= 5     # Number of emission points to consider
#-----------------------------------------------------------------------

# The following are parameters for the perturbation model
h0a  = tpfl[  0  ;  4  ;  8  ;  12  ;  16  ]
σa   = tpfl[  2  ; 1.5 ; 1.8 ;  1.7 ;  1.5 ]
h0i  = tpfl[  150  ;  200  ;  250  ;  300  ;  350  ]
σi   = tpfl[  21   ;  15   ;  18   ;  21   ;  10   ]
	
pfx  = "td"         # Prefix for filename
	
#-----------------------------------------------------------------------

# Fractional uncertainty of atmospheric n due to temperature and 
# pressure uncertainties
δ1	= 0.001

# Fractional uncertainty of ionospheric n 
δ2	= 0.10

# Atmospheric and ionospheric perturbation model
Patm 	= h->metric.P(h,h0a,σa)
Pion 	= h->metric.P(h,h0i,σi)

# Perturbed metric
gp 	    = x->metric.gp(x,δ1,δ2,Patm,Pion)

# Run evaluation function
tdL = squirrel.seval.main(tc,squirrel.locator,gp,Neval,tpfl,tol,ξ,nb,ne)
	
sfx	= "n"*string(ne)*"p"*string(Int(round(δ2*100)))    # Filename suffix

tdLtup = squirrel.seval.td2tup( tdL )
	
Serialization.serialize(dir*pfx*"-"*Nes*"-"*sfx*sufx,tdLtup)
	
#-----------------------------------------------------------------------
	
δ1	= 0.001
δ2	= 0.01
gp 	= x->metric.gp(x,δ1,δ2)
tdS = squirrel.seval.main(tc,squirrel.locator,gp,Neval,tpfl,tol,ξ,nb,ne)
	
sfx	= "n"*string(ne)*"p"*string(Int(round(δ2*100)))

tdStup = squirrel.seval.td2tup( tdS )
	
Serialization.serialize(dir*pfx*"-"*Nes*"-"*sfx*sufx,tdStup)
	
#-----------------------------------------------------------------------
	
td0 = squirrel.seval.main(tc,squirrel.locator,g,Neval,tpfl,tol,ξ,nb,ne)
	
sfx	= "n"*string(ne)*"p0"

td0tup = squirrel.seval.td2tup( td0 )
	
Serialization.serialize(dir*pfx*"-"*Nes*"-"*sfx*sufx,td0tup)
	
#-----------------------------------------------------------------------

tdk = squirrel.seval.main(tck,squirrel.locator,gk,Neval,tpfl,tol,ξ,nb,ne)
	
sfx	= "n"*string(ne)*"k"

tdktup = squirrel.seval.td2tup( tdk )

Serialization.serialize(dir*pfx*"-"*Nes*"-"*sfx*sufx,tdktup)
