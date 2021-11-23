#-----------------------------------------------------------------------
#		EVALUATION SCRIPT
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
tol	= 1e-10
ξ1	= 1e-20
ξ2	= 1e1	
sufx 	= ""
dir 	= "../res/"
Ns	= string(N)
	
#-----------------------------------------------------------------------
#		DESERIALIZE GENERATED SAMPLE FILES (ANALOGUE GEOMETRY)
#-----------------------------------------------------------------------
	
pfx	= "tct" ;
tctloc	= dir*pfx*"-"*Ns*sufx ;

tct	= Serialization.deserialize(tctloc)    ;
	
N    = tct[1] ;
par  = tct[2] ;
X 	 = tct[3] ;
Xtar = tct[4] ;
np   = size(X[1])[2] ;
	
tc	= evalsq.TestCases(par,N,np,X,Xtar) ;

#-----------------------------------------------------------------------
#		DESERIALIZE GENERATED SAMPLE FILES (KERR GEOMETRY)
#-----------------------------------------------------------------------
	
pfx	= "tck" ;
tckloc	= dir*pfx*"-"*Ns*sufx ;

tctk	= Serialization.deserialize(tckloc) ;
	
N	 = tctk[1] ;
par  = tctk[2] ;
X 	 = tctk[3] ;
Xtar = tctk[4] ;
np 	 = size(X[1])[2] ;
	
tck	 = evalsq.TestCases(par,N,np,X,Xtar) ;

#-----------------------------------------------------------------------
ne	= 6     # Number of emission points to consider
#-----------------------------------------------------------------------
	
tpfl = Float64      # Floating point type (determines machine precision)

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
tdL = evalsq.main(tc,squirrel.mlocator,gp,nb,tol,-1,tpfl,true,ξ1,ξ2,ne)
	
sfx	= "n"*string(ne)*"p"*string(Int(round(δ2*100)))    # Filename suffix
	
Serialization.serialize(dir*pfx*"-"*Ns*"-"*sfx*sufx,tdL) # Write to file
	
#-----------------------------------------------------------------------
	
δ1	= 0.001
δ2	= 0.01
gp 	= x->metric.gp(x,δ1,δ2)
tdS = evalsq.main(tc,squirrel.mlocator,gp,nb,tol,-1,Float64,true,ξ1,ξ2,
                  ne)
	
sfx	= "n"*string(ne)*"p"*string(Int(round(δ2*100)))
	
Serialization.serialize(dir*pfx*"-"*Ns*"-"*sfx*sufx,tdS)
	
#-----------------------------------------------------------------------
	
td0 = evalsq.main(tc,squirrel.mlocator,g,nb,tol,-1,Float64,true,ξ1,ξ2,
                  ne)
	
sfx	= "n"*string(ne)*"p0"
	
Serialization.serialize(dir*pfx*"-"*Ns*"-"*sfx*sufx,td0)
	
#-----------------------------------------------------------------------
	
tdk = evalsq.main(tck,squirrel.mlocator,gk,nb,tol,-1,Float64,true,ξ1,ξ2,
                  ne)
	
sfx	= "n"*string(ne)*"k"
	
Serialization.serialize(dir*pfx*"-"*Ns*"-"*sfx*sufx,tdk)
	
#-----------------------------------------------------------------------
	
td6	= (tdL,tdS,td0,tdk)
	
#-----------------------------------------------------------------------
ne	= 5     # Number of emission points to consider
#-----------------------------------------------------------------------
	
tpfl =Float64
	
h0a  = tpfl[  0  ;  4  ;  8  ;  12  ;  16  ]
σa   = tpfl[  2  ; 1.5 ; 1.8 ;  1.7 ;  1.5 ]
h0i  = tpfl[  150  ;  200  ;  250  ;  300  ;  350  ]
σi   = tpfl[  21   ;  15   ;  18   ;  21   ;  10   ]
	
pfx  = "td"
	
#-----------------------------------------------------------------------
	
δ1	= 0.001
δ2	= 0.10
Patm 	= h->metric.P(h,h0a,σa)
Pion 	= h->metric.P(h,h0i,σi)
gp 	    = x->metric.gp(x,δ1,δ2,Patm,Pion)
tdL 	= evalsq.main(tc,squirrel.mlocator,gp,nb,tol,-1,Float64,true,ξ1,
                      ξ2,ne)
	
sfx	= "n"*string(ne)*"p"*string(Int(round(δ2*100)))
	
Serialization.serialize(dir*pfx*"-"*Ns*"-"*sfx*sufx,tdL)
	
#-----------------------------------------------------------------------
	
δ1	= 0.001
δ2	= 0.01
gp 	= x->metric.gp(x,δ1,δ2)
tdS 	= evalsq.main(tc,squirrel.mlocator,gp,nb,tol,-1,Float64,true,ξ1,
                      ξ2,ne)
	
sfx	= "n"*string(ne)*"p"*string(Int(round(δ2*100)))
	
Serialization.serialize(dir*pfx*"-"*Ns*"-"*sfx*sufx,tdS)
	
#-----------------------------------------------------------------------
	
td0 	= evalsq.main(tc,squirrel.mlocator,g,nb,tol,-1,Float64,true,ξ1,
                      ξ2,ne)
	
sfx	= "n"*string(ne)*"p0"
	
Serialization.serialize(dir*pfx*"-"*Ns*"-"*sfx*sufx,td0)
	
#-----------------------------------------------------------------------
	
tdk = evalsq.main(tck,squirrel.mlocator,gk,nb,tol,-1,Float64,true,ξ1,
                  ξ2,ne)
	
sfx	= "n"*string(ne)*"k"
	
Serialization.serialize(dir*pfx*"-"*Ns*"-"*sfx*sufx,tdk)
	
#-----------------------------------------------------------------------
	
td5	= (tdL,tdS,td0,tdk)
	
#-----------------------------------------------------------------------
ne	= 4
#-----------------------------------------------------------------------
	
tpfl =Float64
	
h0a  = tpfl[  0  ;  4  ;  8  ;  12  ;  16  ]
σa   = tpfl[  2  ; 1.5 ; 1.8 ;  1.7 ;  1.5 ]
h0i  = tpfl[  150  ;  200  ;  250  ;  300  ;  350  ]
σi   = tpfl[  21   ;  15   ;  18   ;  21   ;  10   ]
	
pfx  = "td"
	
#-----------------------------------------------------------------------
	
δ1	= 0.001
δ2	= 0.10
Patm 	= h->metric.P(h,h0a,σa)
Pion 	= h->metric.P(h,h0i,σi)
gp 	= x->metric.gp(x,δ1,δ2,Patm,Pion)
tdL 	= evalsq.main(tc,squirrel.slocator,gp,nb,tol,-1,Float64,false,
                      ξ1,ξ2,ne)
	
sfx	= "n"*string(ne)*"p"*string(Int(round(δ2*100)))
	
Serialization.serialize(dir*pfx*"-"*Ns*"-"*sfx*sufx,tdL)
	
#-----------------------------------------------------------------------
	
δ1	= 0.001
δ2	= 0.01
gp 	= x->metric.gp(x,δ1,δ2)
tdS 	= evalsq.main(tc,squirrel.slocator,gp,nb,tol,-1,Float64,false,
                      ξ1,ξ2,ne)
	
sfx	= "n"*string(ne)*"p"*string(Int(round(δ2*100)))
	
Serialization.serialize(dir*pfx*"-"*Ns*"-"*sfx*sufx,tdS)
	
#-----------------------------------------------------------------------
	
td0 	= evalsq.main(tc,squirrel.slocator,g,nb,tol,-1,Float64,false,
                      ξ1,ξ2,ne)
	
sfx	= "n"*string(ne)*"p0"
	
Serialization.serialize(dir*pfx*"-"*Ns*"-"*sfx*sufx,td0)
	
#-----------------------------------------------------------------------
	
tdk = evalsq.main(tck,squirrel.slocator,gk,nb,tol,-1,Float64,false,ξ1,
                  ξ2,ne)
	
sfx	= "n"*string(ne)*"k"
	
Serialization.serialize(dir*pfx*"-"*Ns*"-"*sfx*sufx,tdk)
	
#-----------------------------------------------------------------------
	
td4	= (tdL,tdS,td0,tdk)
	
#-----------------------------------------------------------------------
