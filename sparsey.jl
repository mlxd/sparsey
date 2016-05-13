module sparsey

export idx
export V

# Creates the spatial grid given maximum values and number of samples
function SpatialGrid()

end

# Generate the potential operators along each respective dimension
function V(Grid, Ng::Array{UInt64,1}, Omega::Array{Float64,1}, m::Float64)
	dim = length(Ng)
	V = zeros(dim,)
	Vsp = spzeros(dim,max(Ng))
	for ii = 1:dim
		for jj = Ng[ii]
			Vsp(jj) = 0.5*m*Grid[ii][jj]^2*Omega[ii]^2;
		end
	end
	return Vsp
end
 
#Determines the 1d indexing to an N-d array
function idx(d::Array{Int64},ind::Array{Int64})
	val = 0
	N=length(d)
	for ii = collect(1:N)
		val += (ind[ii]-1)
		if ii == N
			break
		end
		val *= d[ii+1]
	end
	return val+1
end

end
