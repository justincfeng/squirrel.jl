using LinearAlgebra

const RealVec{T<:Real} = Array{T,1}

function zrot( v::RealVec )
	if length(v)==3
		rv=sqrt(v[1]^2+v[2]^2+v[3]^2)		
		Rv=sqrt(v[1]^2+v[2]^2)
		Cosϑ=v[3]/rv
		Sinϑ=Rv/rv
		Cosφ=v[1]/Rv
		Sinφ=v[2]/Rv
		return [ Cosϑ*Cosφ Cosϑ*Sinφ -Sinϑ ; -Sinφ Cosφ 0 ; Sinϑ*Cosφ Sinϑ*Sinφ Cosϑ ]
	else
		print("Not a 3D vector.")	
	end
end	# End zrot
