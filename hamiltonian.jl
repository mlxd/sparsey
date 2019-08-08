module Sparsey

export SpatialGrid
export VHarm
export kronSum
export l2d

import SparseArrays
import Arpack
import LinearAlgebra
import Plots

# Create diagonal using individual arrays from grid. 
# Reshape first dimension diagonal from output into N dimensional grid using Ng values if needed
function kronSum(Ng::Array{Int64,1}, Grid::Array{Array{Float64}})
    dim = length(Ng)
    spd = Array{SparseMatrixCSC{Float64,Int64}}(undef,dim)
    spe = Array{SparseMatrixCSC{Float64,Int64}}(undef,dim)
    for ii=1:(dim)
        spd[ii] = spdiagm(0=>Grid[ii])
        spe[ii] = SparseMatrixCSC{Float64}(I, Ng[ii], Ng[ii])
    end
    for ii=2:dim
        spd[1] = kron(spd[1],spe[ii]) + kron(spe[1],spd[ii])
        #spe[1] = SparseMatrixCSC{Float64}(I, size(spd[1])[1], size(spd[1])[1])
    end
    return spd[1]
end

# 2D Stencil for k grid
function l2d(Ng::Array{Int64,1})
    dim = length(Ng)
    sp = SparseMatrixCSC{Float64}(I, Ng[1], Ng[1])
    E = sparse(2:Ng[1], 1:Ng[1]-1, 1, Ng[1], Ng[1])
    D = E+E' -2*sp;
    K = kron(D,sp) + kron(sp,D)
    return K
end

# Creates the spatial grid given maximum values and number of samples
function SpatialGrid(Ng::Array{Int64,1},GridMax::Array{Float64,1},GridMin::Array{Float64,1})
    dim = length(Ng)
    grid = Array{Array{Float64}}(undef,dim)
    for ii = 1:dim
        grid[ii] = collect(range(GridMin[ii], stop=GridMax[ii], length=Ng[ii]))
    end
    return grid
end

# Generate the potential operators along each respective dimension
function VHarm(Grid::Array{Array{Float64}}, Ng::Array{Int64,1}, Omega::Array{Float64,1}, m::Float64)
    dim = length(Ng)
    V = Array{Array{Float64}}(undef,dim)
    for ii = 1:dim
        V[ii] = zeros(Ng[ii]);
        for jj = 1:Ng[ii]
            V[ii][jj] = 0.5*m*Grid[ii][jj]^2*Omega[ii]^2;
        end
    end
    return V
end

end
