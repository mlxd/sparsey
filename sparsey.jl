module sparsey

export idx
export SpatialGrid
export VHarm
export indicesND
export densify
export kronSum
export derivFD
export l2d

function diagNabbit(M::SparseMatrixCSC{Float64,Int64})
	eigs(M,nev=10,which="SM")
end

# Create diagonal using individual arrays from grid. Reshape first dimension diagonal from output into N dimensional grid using Ng values if needed
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
	return spd[1]
end

function l2d(Ng::Array{Int64,1})
	dim = length(Ng)
	sp = speye(Ng[1])
	E = sparse(2:Ng[1], 1:Ng[1]-1,1,Ng[1],Ng[1])
	D = E+E' -2*sp;
	K = kron(D,sp) + kron(sp,D)
	return K
end

# Generate matrix for n-th derivative of m-th variable. To be finished.
function derivFD(Ng::Array{Int64,1},Grid::Array{Float64},order::Int64)
	dim = length(Ng)
	D = []
	for ii=1:dim
		for jj=1:order
			I = speye(Ng[ii])
			stencil = sparse(2:Ng[ii],1:Ng[ii]-1,1,Ng[ii],Ng[ii])
			stencil += transpose(stencil)
			stencil -= 2*I
			println(stencil)
			if(jj==1)
				println(kron(stencil,I))
				println(kron(I,stencil))
				D = kron(I,stencil)
			elseif(jj==2)
				D += kron(stencil,I)
			else
				D = kron(D,I) + kron(I,D)
			end
		end
	end
	return D
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

#= Pass size of dimensions, and linear index, return i,j,k... value
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
=#

#= Determines the 1d indexing to an N-d array
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
=#

# Return prod(Ng) matrix of vectors in Grid
function densify(Ng::Array{Int64,1},Grid::Array{Array{Float64}})
	idxMax = prod(Ng)
	dim = length(Ng)

end

# Return prod(Ng) matrix from diagonal values in sparse spd
function densify(Ng::Array{Int64,1},spd::SparseMatrixCSC{Float64,Int64})
	idxMax = prod(Ng)
	dim = length(Ng)

end

end
