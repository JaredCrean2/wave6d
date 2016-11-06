# wave equation solver in up to 6 dimensions

module WaveSolver

export runcase

using MPI
include("setup.jl")
include("rk4.jl")
include("communication.jl")

# include generated code
include("generated/loops_5.jl")
include("generated/kernel_5.jl")
include("generated/buffers_2.jl")
include("generated/ic.jl")

function runcase(fname)
  
  f = open(fname, "r")
  nlines  = countlines(f)
  @assert (nlines-1) % 3 == 0
  ndims = div(nlines, 3)
  close(f)

  # countlines leaves the iterator at the end of the file, so reset it
  f = open(fname, "r")

  # the format of the input file must be:
  #  the first nlines lines must contain integers specifying the number of 
  #  grid points in each direction
  #  the first 2*nlines must be floating point numbers specifying the minimum
  #  and maximum coordinate values for each dimension (one per line, alternating
  #  minimum and maximum
  #  The final line must be the a floating point number for the maximum time

  Ns_global = Array(Int, ndims)
  xLs = Array(Float64, 2, ndims)

  println("getting dimension")
  for i=1:ndims
    println("i = ", i)
    str = readline(f)
    println("str = ", str)
    Ns_global[i] = parse(Int, str)
  end

  for i=1:ndims
    str = readline(f)
    xLs[1, i] = parse(Float64, str)
    str = readline(f)
    xLs[2, i] = parse(Float64, str)
  end

  str = readline(f)
  tmax = parse(Float64, str)


  nghost = 2
  params = ParamType(Ns_global, xLs, nghost)

  size_bytes = prod(params.Ns_total_local)*sizeof(Float64)
  println("size of array = ", size_bytes/(1024*1024), " Mbytes")

  dims = zeros(Int, ndims + 1)
  dims[1:end-1] = params.Ns_total_local
  dims[end] = 2
  u_i = Array(Float64, params.Ns_total_local...)

  # apply IC
  IC1(params, u_i)

  # timestep
  tfinal = rk4(step, tmax, u_i, params)

  max_err = calcErr1(params, u_i, tfinal)
  println("max_err = ", max_err)

  MPI.Finalize()


end

function step(params::ParamType, u_i, u_ip1, t)
# single timestep
  startComm(params, u_i)
  finishComm(params, u_i)

  simpleLoop5(params, u_i, u_ip1)

  return nothing
end



end  # module
