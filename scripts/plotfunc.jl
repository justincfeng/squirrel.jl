#-----------------------------------------------------------------------
"""
        rms( x::Vector )

Computes the `rms` value for a vectors.
"""
rms=x->sqrt(dot(x,x)/length(x))            # Root Mean Square calculator

#-----------------------------------------------------------------------
"""
        PrintErr(  h::RealVec , v::RealVec  , e::RealVec , thresh::Tuple , pf::String , units::String )

This function prints the errors, their rms and 95% CL values, and the 
maximum error.

"""
function PrintErr(  h::RealVec , v::RealVec  , e::RealVec , 
                    thresh::Tuple , pf::String , units::String )
    hrms = rms(h)
    vrms = rms(v)
    erms = rms(e)

    n=length(h)

    h95 = h[Int(0.05*n)]
    v95 = v[Int(0.05*n)]
    e95 = e[Int(0.05*n)]

    hcthr = count(x->x>thresh[1], h)
    vcthr = count(x->x>thresh[2], v)
    ecthr = count(x->x>thresh[3], e)

    open(pf*"-res.txt","a") do io
        println(io,"hrms =","\t",hrms,"\t",
                "# Horizontal error rms (",units,")")
        println(io,"vrms =","\t",vrms,"\t",
                "# Vertical error rms (",units,")")
        println(io,"erms =","\t",erms,"\t",
                "# Full error rms (",units,")")
        println(io,"h95 =","\t",h95,"\t",
                "# Horizontal error 95% CL (",units,")")
        println(io,"v95 =","\t",v95,"\t",
                "# Vertical error 95% CL (",units,")")
        println(io,"e95 =","\t",e95,"\t",
                "# Full error 95% CL (",units,")")
        println(io,"hcthr =","\t",hcthr,"\t",
                "# Horizontal errors above threshold of ",thresh[1]," ",
                units)
        println(io,"vcthr =","\t",vcthr,"\t",
                "# Vertical errors above threshold of ",thresh[2]," ",
                units)
        println(io,"ecthr =","\t",ecthr,"\t",
                "# Full errors above threshold of ",thresh[3]," ",units)
        println(io,"hmax =","\t",h[1],"\t",
                "# Largest horizontal error (",units,")")
        println(io,"vmax =","\t",v[1],"\t",
                "# Largest vertical error (",units,")")
        println(io,"emax =","\t",e[1],"\t",
                "# Largest total error (",units,")")
    end

    return (hrms,vrms,erms,h95,v95,e95,hcthr,vcthr,ecthr)
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
PlotErr(   err::RealVec  , lines::Tuple , histpar::Tuple , xr::Tuple , yr::Tuple , labels::Tuple , plotname::String )

This function plots the errors.

"""
function PlotErr(   err::RealVec  , lines::Tuple , histpar::Tuple , 
                    xr::Tuple , yr::Tuple , labels::Tuple , 
                    plotname::String )
    plot!(vline([lines[1],lines[2]]))
    histogram!(
        err ,
        alpha=0.50 ,
        bins=histpar[2]:histpar[1]:histpar[3],
        xlims=xr ,
        ylims=yr ,
        xlabel=labels[1] ,
        ylabel=labels[2] ,
        leg=false ,
        fontfamily="Times"
    )
    savefig(plotname)
end     #---------------------------------------------------------------
