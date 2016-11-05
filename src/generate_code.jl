# generate all the files

include("generate_kernel.jl")
include("generate_loops.jl")
include("generate_buffer.jl")

# input values
stencil = ["1.0", "2.0", "3.0", "4.0", "5.0"]
maxdim = 5
neq = 2  # number of equations
prefix = "generated/"


npts = length(stencil)
nghost = div(npts, 2)

generate_kernel(maxdim, stencil, neq, prefix)
generate_loops(maxdim, npts, prefix)
generate_buffers(2, maxdim, prefix)
