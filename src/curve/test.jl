include("hilbert.jl")
#=
ndims = 6
npoints = 16^6

checkDimensions(ndims, 1024)

coords, idxs = getState(ndims, npoints)

@time loadNpoints(npoints, ndims, coords, idxs)
@time loadNpoints(npoints, ndims, coords, idxs)

writedlm("hdata3.dat", idxs.' - 1)

nblocks = getNumBlocks(1024, 2, 64)
println("nblocks = ", nblocks)

ndims = 3
coords, idxs = getState(ndims, npoints)
checkDimensions(ndims, 16)

loadNpoints(npoints, ndims, coords, idxs)

writedlm("hdata3.dat", idxs.' - 1)
=#

ndims = 2
npoints = 8^ndims
coords, idxs = getState(ndims, npoints)
loadNpoints(npoints, ndims, coords, idxs)
writedlm("hdata2.dat", idxs.')

ndims = 3
npoints = 8^ndims
coords, idxs = getState(ndims, npoints)
loadNpoints(npoints, ndims, coords, idxs)
writedlm("hdata3.dat", idxs.')

