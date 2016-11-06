# generate all the files

include("generate_kernel.jl")
include("generate_loops.jl")
include("generate_buffer.jl")
include("generate_ic.jl")
include("generate_calcerr.jl")

# input values
stencil = ["-1/12", "4/3", "-5/2", "4/3", "-1/12"]
maxdim = 6
neq = 2  # number of equations
prefix = "generated/"


npts = length(stencil)
nghost = div(npts, 2)

generate_kernel(maxdim, stencil, neq, prefix)
generate_loops(maxdim, npts, prefix)
generate_buffers(2, maxdim, prefix)
generate_ic(maxdim, prefix)
generate_calcerr(maxdim, prefix)
