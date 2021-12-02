#---------------------------------------------------------------------------------------;
#		PLOTS AND RESULTS
#---------------------------------------------------------------------------------------;

using Plots, LinearAlgebra, Serialization, BenchmarkTools
using CoordinateTransformations, Geodesy

include("../src/type.jl")
include("../src/evalsq.jl")
include("../src/srl5.jl")
include("../src/squirrel.jl")
include("../src/geocoord.jl")
include("../src/metric.jl")
include("plotfunc.jl")

g  	= metric.g
gk 	= metric.ge

N	= 100000
# N	= 100000

nb	= 24
tol	= 1e-10
ξ1	= 1e-18
ξ2	= 1e1	
sufx 	= ""
dir 	= "../res/"
Ns	= string(N)

rms=x->sqrt(dot(x,x)/length(x))
mfact   = 0.004435
cmfact  = 0.4435
mmfact  = 4.435

#---------------------------------------------------------------------------------------;
ne	= 6
#---------------------------------------------------------------------------------------;
	
tpfl =Float64
	
h0a  = tpfl[  0  ;  4  ;  8  ;  12  ;  16  ]
σa   = tpfl[  2  ; 1.5 ; 1.8 ;  1.7 ;  1.5 ]
h0i  = tpfl[  150  ;  200  ;  250  ;  300  ;  350  ]
σi   = tpfl[  21   ;  15   ;  18   ;  21   ;  10   ]

Patm 	= h->metric.P(h,h0a,σa)
Pion 	= h->metric.P(h,h0i,σi)
	
pfx  = "td"
	
#---------------------------------------------------------------------------------------;

δ1	= 0.001
δ2	= 0.10
gp 	= x->metric.gp(x,δ1,δ2,Patm,Pion)
sfx	= "n"*string(ne)*"p"*string(Int(round(δ2*100)))

tdL 	= Serialization.deserialize(dir*pfx*"-"*Ns*"-"*sfx*sufx)

lfact   = mfact

thresh = (20,20,20)

h   = lfact*sort(tdL.erh,rev=true)
v   = lfact*sort(tdL.erv,rev=true)
e   = lfact*sort(tdL.err,rev=true)

pf = dir*pfx*"-"*Ns*"-"*sfx*sufx

ers = PrintErr( h , v , e , thresh , pf , "m" )

histpar  = (0.05,0.0,50)
xrange   = (0,10)
yrange   = (0,15000)
labelsh  = ("Horizontal error (m)","Number of test cases")
labelsv  = ("Vertical error (m)","Number of test cases")

PlotErr( h , (ers[1],ers[4]) , histpar , xrange , yrange , labelsh 
        , dir*"plotHerr-"*sfx*".pdf" )
PlotErr( v , (ers[2],ers[5]) , histpar , xrange , yrange , labelsv 
        , dir*"plotVerr-"*sfx*".pdf" )

#---------------------------------------------------------------------------------------;

δ1	= 0.001
δ2	= 0.01
gp 	= x->metric.gp(x,δ1,δ2)
sfx	= "n"*string(ne)*"p"*string(Int(round(δ2*100)))

tdS 	= Serialization.deserialize(dir*pfx*"-"*Ns*"-"*sfx*sufx)

lfact   = cmfact

thresh = (200,200,200)

h   = lfact*sort(tdS.erh,rev=true)
v   = lfact*sort(tdS.erv,rev=true)
e   = lfact*sort(tdS.err,rev=true)

pf = dir*pfx*"-"*Ns*"-"*sfx*sufx

ers = PrintErr( h , v , e , thresh , pf , "cm" )

histpar  = (0.5,0.0,500)
xrange   = (0,100)
yrange   = (0,15000)
labelsh  = ("Horizontal error (cm)","Number of test cases")
labelsv  = ("Vertical error (cm)","Number of test cases")

PlotErr( h , (ers[1],ers[4]) , histpar , xrange , yrange , labelsh 
        , dir*"plotHerr-"*sfx*".pdf" )
PlotErr( v , (ers[2],ers[5]) , histpar , xrange , yrange , labelsv 
        , dir*"plotVerr-"*sfx*".pdf" )

#---------------------------------------------------------------------------------------;

sfx	= "n"*string(ne)*"p0"

td0 	= Serialization.deserialize(dir*pfx*"-"*Ns*"-"*sfx*sufx)

lfact   = mmfact

thresh = (30,30,30)

h   = lfact*sort(td0.erh,rev=true)
v   = lfact*sort(td0.erv,rev=true)
e   = lfact*sort(td0.err,rev=true)

pf = dir*pfx*"-"*Ns*"-"*sfx*sufx

ers = PrintErr( h , v , e , thresh , pf , "mm" )

histpar  = (0.01,0.0,1000)
xrange   = (0,2)
yrange   = (0,20000)
labelsh  = ("Horizontal error (cm)","Number of test cases")
labelsv  = ("Vertical error (cm)","Number of test cases")

PlotErr( h , (ers[1],ers[4]) , histpar , xrange , yrange , labelsh 
        , dir*"plotHerr-"*sfx*".pdf" )
PlotErr( v , (ers[2],ers[5]) , histpar , xrange , yrange , labelsv 
        , dir*"plotVerr-"*sfx*".pdf" )


#---------------------------------------------------------------------------------------;

pfx  = "tck"

sfx	= "n"*string(ne)*"k"

tdk = Serialization.deserialize(dir*pfx*"-"*Ns*"-"*sfx*sufx)

lfact   = mmfact

thresh = (20,20,20)

h   = lfact*sort(tdk.erh,rev=true)
v   = lfact*sort(tdk.erv,rev=true)
e   = lfact*sort(tdk.err,rev=true)

hC  = cmfact*sort(tdk.erhC,rev=true)
vC  = cmfact*sort(tdk.ervC,rev=true)
eC  = cmfact*sort(tdk.errC,rev=true)

pf = dir*pfx*"-"*Ns*"-"*sfx*sufx

ers = PrintErr( h , v , e , thresh , pf , "mm" )
ersC = PrintErr( hC , vC , eC , thresh , pf*"C" , "cm" )

histpar  = (0.0025,0.0,100)
xrange   = (0,0.5)
yrange   = (0,30000)
labelsh  = ("Horizontal error (mm)","Number of test cases")
labelsv  = ("Vertical error (mm)","Number of test cases")

PlotErr( h , (ers[1],ers[4]) , histpar , xrange , yrange , labelsh 
        , dir*"plotHerr-"*sfx*".pdf" )
PlotErr( v , (ers[2],ers[5]) , histpar , xrange , yrange , labelsv 
        , dir*"plotVerr-"*sfx*".pdf" )

histpar  = (0.0025,0.0,100)
xrange   = (0,0.5)
yrange   = (0,30000)
labelsh  = ("Horizontal error (cm)","Number of test cases")
labelsv  = ("Vertical error (cm)","Number of test cases")

PlotErr( h , (ersC[1],ersC[4]) , histpar , xrange , yrange , labelsh 
        , dir*"plotHerr-"*sfx*"C.pdf" )
PlotErr( v , (ersC[2],ersC[5]) , histpar , xrange , yrange , labelsv 
        , dir*"plotVerr-"*sfx*"C.pdf" )

#---------------------------------------------------------------------------------------;


