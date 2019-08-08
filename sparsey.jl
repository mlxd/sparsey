module Sparsey
include("hamiltonian.jl")
export diag

struct HOscParams
    omegaXY::Array{Float64}
    mass::Float64
    gridDivs::Array{Integer}
    gridScale::Array{Float64}
end

function diag(osc::HOscParams, numEig::Integer)
    # Set number of eigenvalues, and spatial dim of grid

    Ngx = [osc.gridDivs[1], osc.gridDivs[2]]
    dim = length(Ngx)
    gMax = [osc.gridScale[1], osc.gridScale[2]]
    gMin = [-osc.gridScale[1], -osc.gridScale[2]]

    g0 = SpatialGrid(Ngx, gMax, gMin)
    V = VHarm(g0, Ngx, [osc.omegaXY[1], osc.omegaXY[2]], osc.mass)

    # Returns sparse diagonal matrix with linear indexed (lexicographical) 
    # values of the vectors V, with dimensions Ngx
    VD=kronSum(Ngx,V) 

    dx = abs.(g0[1][1] - g0[1][2])
    KD = -(1./(2.*osc.mass))*l2d(Ngx)./dx^2

    H = VD + KD
    val,vec = eigs(H, nev=numEig, which=:SM); #Calculate 8 of the smallest eigenvectors
    return val, vec
end


end
