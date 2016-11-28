# generate all the files

include("generate_kernel.jl")
include("generate_loops.jl")
include("generate_loops_blocked.jl")
include("generate_loops_hilbert.jl")
include("generate_loops_blockhilbert.jl")
include("generate_buffer.jl")
include("generate_ic.jl")
include("generate_calcerr.jl")
include("generate_step.jl")

# input values
stencil = ["1/12", "-2/3", "0", "2/3", "-1/12"]
maxdim = 6
neq = 2  # number of equations
prefix = "generated/"
#blocksizes = [2, 4, 8, 16]
blocksizes = [2, 4, 8, 16]

# select whether or not to tell the compiler that vectorization is legal
vectorize = false
if vectorize
  macro_name = "@simd"
else
  macro_name = ""
end

npts = length(stencil)
nghost = div(npts, 2)

generate_kernel(maxdim, stencil, neq, prefix)
generate_loops(maxdim, npts, prefix)
for bs in blocksizes
  generate_loops_blocked(maxdim, npts, bs, prefix)
  generate_loops_blockhilbert(maxdim, npts, bs, prefix)
end
generate_loops_hilbert(maxdim, npts, prefix)
generate_buffers(2, maxdim, prefix)
generate_ic(maxdim, prefix)
generate_calcerr(maxdim, prefix)
generate_step(maxdim, npts, blocksizes, prefix)

