module sparsey

export idx
export SpatialGrid
export VHarm
export indicesND
export densify
export kronSum

# Create diagonal using individual arrays from grid
function kronSum(Ng::Array{Int64,1},Grid::Array{Array{Float64}})
	dim = length(Ng)
	spd = Array{SparseMatrixCSC{Float64,Int64}}(dim)
	spe = Array{SparseMatrixCSC{Float64,Int64}}(dim)
	for ii=1:(dim)
		spd[ii] = spdiagm(Grid[ii],0,Ng[ii],Ng[ii])
		spe[ii] = speye(Ng[ii])
	end
	for ii=2:dim
		spd[1] = kron(spd[1],spe[ii]) + kron(spe[1],spd[ii])
		spe[1] = speye(size(spd[1])[1])
	end
	return spd
end

# Creates the spatial grid given maximum values and number of samples
function SpatialGrid(Ng::Array{Int64,1},GridMax::Array{Float64,1},GridMin::Array{Float64,1})
	dim = length(Ng)
	grid = Array{Array{Float64}}(dim)
	for ii = 1:dim
		grid[ii] = collect(linspace(GridMin[ii],GridMax[ii],Ng[ii]))
	end
	return grid
end

# Generate the potential operators along each respective dimension
function VHarm(Grid::Array{Array{Float64}}, Ng::Array{Int64,1}, Omega::Array{Float64,1}, m::Float64)
	dim = length(Ng)
	V = Array{Array{Float64}}(dim)
	for ii = 1:dim
		V[ii] = zeros(Ng[ii]);
		for jj = 1:Ng[ii]
			V[ii][jj] = 0.5*m*Grid[ii][jj]^2*Omega[ii]^2;
		end
	end
	return V
end

function densify(Ng::Array{Int64,1},Grid::Array{Array{Float64}})
	idxMax = prod(Ng)
	dim = length(Ng)
	S = Array{Float64}(idxMax)
	k = Grid[1]
	for ii = 2:dim
		k = kron(k,Grid[ii]) #Double check this. Not sure if it will work here
	end
	println(k)
	return k
end

# Pass size of dimensions, and linear index, return i,j,k... value
function indicesND(Ng::Array{Int64,1},linIdx::Int64)
	dim = length(Ng)
	slices = Array{Int64}(dim)
	for ii = 1:dim
		slices[ii] = prod(Ng[1:dim-ii+1])
	end
	sl = (floor(Integer,linIdx./slices))

	idx = find(x-> x != 0,sl)
	println(idx)

	slp = prod(sl[idx])
	println(slp)

	r = linIdx % slp;
	println(r)
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
