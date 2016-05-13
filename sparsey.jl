module sparsey

export idx
export VHarm

# Creates the spatial grid given maximum values and number of samples
function SpatialGrid()

end

# Generate the potential operators along each respective dimension
function VHarm(Grid, Ng::Array{Int64,1}, Omega::Array{Float64,1}, m::Float64)
	dim = length(Ng)
	V = Array{Array{Float64}}(dim)
	for ii = 1:dim
		V[ii] = zeros(Ng[ii]);
		for jj = 1:Ng[ii]
			V[ii][jj] = 0.5*m*Grid[ii][jj]^2*Omega[ii]^2;
		end
	end
	println(V)
	return V
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
