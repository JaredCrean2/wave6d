include("hilbert.jl")

ndims = 2
npoints = 1024

checkDimensions(ndims, 1024)

coords, idxs = getState(ndims, npoints)

@time loadNpoints(npoints, ndims, coords, idxs)
@time loadNpoints(npoints, ndims, coords, idxs)

writedlm("hdata3.dat", idxs.' - 1)

nblocks = getNumBlocks(1024, 2, 64)
println("nblocks = ", nblocks)
