# wave equation solver in up to 6 dimensions

module WaveSolver

export runcase

using MPI
include("setup.jl")
include("rk4.jl")
include("lserk.jl")
include("communication.jl")

# include generated code
include("generated/loops_5.jl")
include("generated/kernel_5.jl")
include("generated/buffers_2.jl")
include("generated/ic.jl")
include("generated/calcerr.jl")

# These would be preprocessor defins controlled by the build system in C
global const FORCE_SYNC = true  # force synchronization at the beignning of 
                                # every stage
global const USE_LOW_STORAGE = true

function runcase(fname)
  
  Ns_global, xLs, tmax, write_conv = parseinput(fname)

  ndims = length(Ns_global)
  nghost = 2
  params = ParamType(Ns_global, xLs, nghost)

  size_bytes = prod(params.Ns_global + 2*nghost)*sizeof(Float64)
  if params.comm_rank == 0
    println("global size of array = ", size_bytes/(1024*1024), " Mbytes")
  end

  dims = zeros(Int, ndims + 1)
  dims[1:end-1] = params.Ns_total_local
  dims[end] = 2
  u_i = Array(Float64, dims...)
  fill!(u_i, 0.0)

  # apply IC
  IC1(params, u_i)

  MPI.Barrier(params.comm)

  # timestep
  if USE_LOW_STORAGE
    tfinal = lserk(step, tmax, u_i, params)
  else
    tfinal = rk4(step, tmax, u_i, params)
  end


  max_err = calcErr1(params, u_i, tfinal)
  max_err = MPI.Allreduce(max_err, MPI.MAX, params.comm)
  if params.comm_rank == 0
    println("max_err = ", max_err)
    if write_conv == 1
      f2 = open("convergence.dat", "a+")
      println(f2, maximum(params.delta_xs), " ", max_err)
      close(f2)
    end
  end

end

function parseinput(fname)
  f = open(fname, "r")
  nlines  = countlines(f)
  @assert (nlines-2) % 3 == 0
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

  for i=1:ndims
    str = readline(f)
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

  str = readline(f)
  write_conv = parse(Int, str)

  close(f)

  return Ns_global, xLs, tmax, write_conv

end



function step{N}(params::ParamType{N}, u_i, u_ip1, t)
# single timestep
  params.t = t

  startComm(params, u_i)

  finishComm(params, u_i)

  if FORCE_SYNC
    MPI.Barrier(params.comm)
  end


  simpleLoop5(params, u_i, u_ip1)

  params.itr += 1

  return nothing
end



end  # module
