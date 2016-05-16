push!(LOAD_PATH,"/home/mlxd/prog/sparsey")
using sparsey
using PyPlot

Ngx = [100,100]
dim = length(Ngx)

gMax = [10.,10.]
gMin = [-10.,-10.]

g0 = SpatialGrid(Ngx,gMax,gMin)
V = VHarm(g0,Ngx,[1.,1.],1.)

#densify([Ngx],V) #Convert sparse matrix of size prod(Ngx)^2 to dense matrix of size prod(Ngx). Alternatively, if given a vector of the values in each dimension, returns similarly N-d matrix/Arrays
#indicesND(Ngx,3) #Given the dimensions Ngx of the system, and a linear index, returns the values of index in each respective dimension

VD=kronSum(Ngx,V) #Returns sparse diagonal matrix with linear indexed (lexicographical) values of the vectors V, with dimensions Ngx

dx = abs(g0[1][1] - g0[1][2])
KD = -0.5*l2d(Ngx)./dx^2

H = VD + KD
val,vec = eigs(H,nev=8,which="SM") #Calculate 8 of the smallest eigenvectors

plot(val)
#pcolor(reshape(vec[:,1],Ngx[1],Ngx[2]))
#colorbar()
show()
