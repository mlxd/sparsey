push!(LOAD_PATH,"/home/mlxd/prog/sparsey")
using sparsey

Ngx = [8,4,6]
dim = 3

# g = Array{Array{Float64}}(dim)
# g[1] = collect(linspace(-10,10,Ngx[1]))
# g[2] = collect(linspace(-4,7,Ngx[2])) 
# g[3] = collect(linspace(-1,1,Ngx[3]))


gMax = [10.,7.,1.]
gMin = [-10.,-4.,-1.]

g = SpatialGrid(Ngx,gMax,gMin)
VHarm(g,Ngx,[1.,1.,1.],1.)
