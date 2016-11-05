# wave equation solver in up to 6 dimensions

module WaveSolver
using MPI
include("setup.jl")

include("generated/loops_5.jl")
include("generated/kernel_5.jl")
include("generated/buffers_2.jl")

# include the kernels here


end
