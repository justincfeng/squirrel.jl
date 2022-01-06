#-----------------------------------------------------------------------
#		PLOTS AND RESULTS
#-----------------------------------------------------------------------

using Plots, LinearAlgebra, Serialization, BenchmarkTools

include("../src/squirrel.jl")
include("../src/metric.jl")
include("../src/type.jl")
include("plotfunc.jl")

#-----------------------------------------------------------------------
#       PLOT PARAMETERS
#-----------------------------------------------------------------------

Neval	= 1000                  # Evaluated test cases (for filenames)
N	= 1000                  # Number of test cases to plot

pfx     = "td"                  # Filename string prefix
sufx 	= ""                    # Filename string suffix
dir 	= "../res/"             # Results directory
Nes	= string(Neval)         # Turn Neval into a string
Ns      = string(N)             # Turn Ns into a string

tpfl    = Float64               # Floating point precision to use

hisfac  = Int(round(N/100))     # Histogram factor

#-----------------------------------------------------------------------
#       CONVERSION FACTORS
#-----------------------------------------------------------------------

mfact   = 0.004435              # conversion factor to meters
cmfact  = 0.4435                # conversion factor to centimeters
mmfact  = 4.435                 # conversion factor to millimeters

#-----------------------------------------------------------------------
ne	= 5
#-----------------------------------------------------------------------
	
#-----------------------------------------------------------------------

δ1	= 0.001                 # Atmospheric perturbation parameter
δ2	= 0.10                  # Ionospheric perturbation parameter
sfx	= "n"*string(ne)*"p"*string(Int(round(δ2*100))) # Filename str.

tdLtup	= Serialization.deserialize(dir*pfx*"-"*Ns*"-"*sfx*sufx)
tdL     = squirrel.seval.tup2td( tdLtup )

lfact   = mfact                         # Choose conversion factors
unitstr = "m"                           # Unit string

thresh = (20,20,20)                             # Threshold tuples 

h   = lfact*sort(tdL.erh,rev=true)              # horizontal errors
v   = lfact*sort(tdL.erv,rev=true)              # vertical errors
e   = lfact*sort(tdL.err,rev=true)              # total error

pf = dir*pfx*"-"*Ns*"-"*sfx*sufx                # Filename string

ers = PrintErr( h , v , e , thresh , pf , unitstr )     # Error string

binwidth = 0.05                                 # Histogram bin width
histL    = 0.0                                  # Histogram lower limit
histU    = 50.0                                 # Histogram upper limit

xrange   = (0,10)               # x range for histogram plot 
yrange   = (0,15*hisfac)        # y range for histogram plot 

histpar  = (binwidth,histL,histU)
labelsh  = ("Horizontal error ("*unitstr*")","Number of test cases")
labelsv  = ("Vertical error ("*unitstr*")","Number of test cases")

PlotErr( h , (ers[1],ers[4]) , histpar , xrange , yrange , labelsh ,
         dir*"plotHerr-"*sfx*".pdf" )
PlotErr( v , (ers[2],ers[5]) , histpar , xrange , yrange , labelsv ,
         dir*"plotVerr-"*sfx*".pdf" )

#-----------------------------------------------------------------------

δ1	= 0.001                 # Atmospheric perturbation parameter
δ2	= 0.01                  # Ionospheric perturbation parameter
sfx	= "n"*string(ne)*"p"*string(Int(round(δ2*100))) # Filename str.

tdStup	= Serialization.deserialize(dir*pfx*"-"*Ns*"-"*sfx*sufx)
tdS     = squirrel.seval.tup2td( tdStup )

lfact   = cmfact                        # Choose conversion factors
unitstr = "cm"                          # Unit string

thresh = (20,20,20)                          # Threshold tuples 

h   = lfact*sort(tdS.erh,rev=true)              # horizontal errors
v   = lfact*sort(tdS.erv,rev=true)              # vertical errors
e   = lfact*sort(tdS.err,rev=true)              # total error

pf = dir*pfx*"-"*Ns*"-"*sfx*sufx                # Filename string

ers = PrintErr( h , v , e , thresh , pf , unitstr )     # Error string

binwidth = 0.5                                  # Histogram bin width
histL    = 0.0                                  # Histogram lower limit
histU    = 500.0                                # Histogram upper limit

xrange   = (0,100)              # x range for histogram plot 
yrange   = (0,15*hisfac)        # y range for histogram plot 

histpar  = (binwidth,histL,histU)
labelsh  = ("Horizontal error ("*unitstr*")","Number of test cases")
labelsv  = ("Vertical error ("*unitstr*")","Number of test cases")

PlotErr( h , (ers[1],ers[4]) , histpar , xrange , yrange , labelsh ,
         dir*"plotHerr-"*sfx*".pdf" )
PlotErr( v , (ers[2],ers[5]) , histpar , xrange , yrange , labelsv ,
         dir*"plotVerr-"*sfx*".pdf" )

#-----------------------------------------------------------------------

sfx	= "n"*string(ne)*"p0"

td0tup 	= Serialization.deserialize(dir*pfx*"-"*Ns*"-"*sfx*sufx)
td0     = squirrel.seval.tup2td( td0tup )

lfact   = mmfact                        # Choose conversion factors
unitstr = "mm"                          # Unit string

thresh = (20,20,20)                             # Threshold tuples 

h   = lfact*sort(td0.erh,rev=true)              # horizontal errors
v   = lfact*sort(td0.erv,rev=true)              # vertical errors
e   = lfact*sort(td0.err,rev=true)              # total error

pf = dir*pfx*"-"*Ns*"-"*sfx*sufx                # Filename string

ers = PrintErr( h , v , e , thresh , pf , unitstr )     # Error string

binwidth = 0.01                                 # Histogram bin width
histL    = 0.0                                  # Histogram lower limit
histU    = 10.0                                 # Histogram upper limit

xrange   = (0,2)                # x range for histogram plot 
yrange   = (0,20*hisfac)        # y range for histogram plot 

histpar  = (binwidth,histL,histU)
labelsh  = ("Horizontal error ("*unitstr*")","Number of test cases")
labelsv  = ("Vertical error ("*unitstr*")","Number of test cases")

PlotErr( h , (ers[1],ers[4]) , histpar , xrange , yrange , labelsh ,
         dir*"plotHerr-"*sfx*".pdf" )
PlotErr( v , (ers[2],ers[5]) , histpar , xrange , yrange , labelsv ,
         dir*"plotVerr-"*sfx*".pdf" )

#-----------------------------------------------------------------------

sfx	= "n"*string(ne)*"k"

tdktup  = Serialization.deserialize(dir*pfx*"-"*Ns*"-"*sfx*sufx)
tdk     = squirrel.seval.tup2td( tdktup )

lfact   = mmfact                        # Choose conversion factors
unitstr = "mm"                          # Unit string

thresh = (20,20,20)                             # Threshold tuples 

h   = lfact*sort(tdk.erh,rev=true)              # horizontal errors
v   = lfact*sort(tdk.erv,rev=true)              # vertical errors
e   = lfact*sort(tdk.err,rev=true)              # total error

pf = dir*pfx*"-"*Ns*"-"*sfx*sufx                # Filename string

ers = PrintErr( h , v , e , thresh , pf , unitstr )     # Error string

binwidth = 0.0025                               # Histogram bin width
histL    = 0.0                                  # Histogram lower limit
histU    = 5.0                                # Histogram upper limit

xrange   = (0,0.5)              # x range for histogram plot 
yrange   = (0,30*hisfac)        # y range for histogram plot  

histpar  = (binwidth,histL,histU)
labelsh  = ("Horizontal error ("*unitstr*")","Number of test cases")
labelsv  = ("Vertical error ("*unitstr*")","Number of test cases")

PlotErr( h , (ers[1],ers[4]) , histpar , xrange , yrange , labelsh 
        , dir*"plotHerr-"*sfx*".pdf" )
PlotErr( v , (ers[2],ers[5]) , histpar , xrange , yrange , labelsv 
        , dir*"plotVerr-"*sfx*".pdf" )

pfC = pf*"C"

lfactC   = cmfact                       # Choose conversion factors
unitstrC = "cm"                         # Unit string

threshC = (2,5,20)                            # Threshold tuples 

hC  = lfactC*sort(tdk.erhC,rev=true)            # horizontal errors
vC  = lfactC*sort(tdk.ervC,rev=true)            # vertical errors
eC  = lfactC*sort(tdk.errC,rev=true)            # total errors

ersC = PrintErr( hC , vC , eC , threshC , pfC , unitstrC )  # Error str

labelsh  = ("Horizontal error ("*unitstrC*")","Number of test cases")
labelsv  = ("Vertical error ("*unitstrC*")","Number of test cases")

binwidth = 0.00025                              # Histogram bin width
histpar  = (binwidth,histL,histU)
xrange   = (0,0.05)             # x range for histogram plot 
yrange   = (0,20*hisfac)        # y range for histogram plot 
PlotErr( hC , (ersC[1],ersC[4]) , histpar , xrange , yrange , labelsh ,
         dir*"plotHerr-"*sfx*"C.pdf" )

binwidth = 0.02                              # Histogram bin width
histpar  = (binwidth,histL,histU)
xrange   = (0,4)               # x range for histogram plot 
yrange   = (0,4*hisfac)        # y range for histogram plot  
PlotErr( vC , (ersC[2],ersC[5]) , histpar , xrange , yrange , labelsv ,
         dir*"plotVerr-"*sfx*"C.pdf" )

#-----------------------------------------------------------------------
