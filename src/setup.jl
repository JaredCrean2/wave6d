type ParamType{N}
  delta_xs::Array{Float64, 1}
  deltax_invs2::Array{Float64, 1}
  comm::MPI.Comm
  Ns_global::Array{Int, 1}  # number of regular grid points in each direction
  Ns_local::Array{Int, 1}
  Ns_total_local::Array{Int, 1}  # number to regular grid points + ghosts in each direction
  ias::Array{Int, 1}
  ibs::Array{Int, 1}
  send_reqs::Array{MPI.Request, 1}
  recv_reqs::Array{MPI.Request, 1}
  peernums::Array{Int, 2}  # 2 x ndim array containing peer numbers
  xLs::Array{Float64, 2}  # xmin and xmax for each dimension
  nghost::Int
  coords::Array{LinSpace{Float64}, 1}
end

include("setup2.jl")  # the auto-generated part

function ParamType(Ns_global::Array{Int, 1}, xLs::Array{Float64, 2}, nghost)
# Ns = number of grid points (not including ghosts
# xls = 2 x ndim array of xmin and xmax for each dimension

  Ns_local = copy(Ns_global)  # change when parallelizing
  Ns_total_local = Ns_local + 2*nghost

  N = length(Ns_global)
  comm = MPI.COMM_WORLD

  ias = nghost*ones(N) + 1
  ibs = zeros(ias)
  for i=1:N
    ibs[i] = Ns_total_local[i] - nghost
  end

  delta_xs = zeros(N)
  for i=1:N
    xmin = xLs[1, i]
    xmax = xLs[2, i]
    N_i = Ns[i]
    delta_xs[i] = (xmax - xmin)/(N_i - 1)
  end

  delta_xinvs2 = 1./(delta_xs.^2)

  send_reqs = Array(MPI.Request, 0)
  recv_reqs = Array(MPI.Request, 0)
  peernums = Array(Int, 0)

  # TODO update this when parallelizing
  coords = Array(LinSpace{Float64}, N)
  for i=1:N
    xmin = xls[1, i]
    xmax = xLs[2, i]
    delta_x = delta_xs[i]
    Ntot = Ns_total_local[i]

    xmin = xmin - nghost*delta_x
    xmax = xmax + nghost*delta_x

    coords[i] = linspace(xmin, xmax, Ntot)
  end

  return ParamType{N}(delta_xs, deltax_invs2, comm, Ns_global, Ns_local, Ns_total_local,
                      ias, ibs, send_reqs, recv_reqs, peernums, xLs, nghost, coords)
end


"""
  This function calculates the number of processes along each axis
  via exhaustive search
"""
function mpiCalculation(comm, N::Integer)
# comm is a MPI communicator, N is the number of dimensions
# computes the number of processes along each axis by solving a 
# least squared problem (brute force)

# it then gets the N indices that describe the location of this rank in
# the grid, and the ranks of its neighbors

  comm_size = MPI.Comm_size(comm)
  comm_rank = MPI.Comm_rank(comm)

  matches, idx = getMpiMatches(comm_size, N)

  # figure out how many rank in each dimension by minimizing
  # the square of the difference from the optimal value.
  # the optimal value is the Nth root of the number of 
  # processes
  max_diff = typemax(Int)
  idx_opt = 0
  optimal_val = Float64(comm_size)^(1/N)

  for i=1:idx
    diff_squared = 0.0
    for j=1:N
      diff = optimal_val - matches[i, j]
      diff_squared += diff*diff
    end

    if diff_squared < max_diff
      max_diff = diff_squared
      idx_opt = i
    end
  end


  dims = (matches[idx_opt, :]...)
  # generate an N dimensional array and assign MPI
  # ranks to it

  # this is type-unstable, but that's fine
  # TODO: see if it is possible to avoid allocating this

  # the calculations done below are equivalent:
  #=
  rankgrid = zeros(Int, dims...)

  for i=0:(comm_size-1)
    rankgrid[i] = i
  end
  =#
  # where rankgrid is the array that maps ranks to their position on the grid

  # get i, j, k... that describe this ranks location in the grid
  my_subs = ind2sub(dims, comm_rank + 1)
  neighbor_subs = copy(my_subs)
  
  # get the rank number of the adjacent blocks for each dimension
  peer_nums = Array(Int, 2, N)

  for i=1:N
    idx_i = my_subs[i]
    # this behaves correctly even if a dimension is 1
    # left rank
    if idx_i == 1  # this is the leftmost process
      left_idx = dims[i]  # the maximum process
    else
      left_idx = idx - 1
    end

    # right rank
    if idx_i == dims[i]
      right_idx = 1
    else
      right_idx = idx_i + 1
    end


    # substitute left_idx, right_idx into my_subs to compute the rank
    copy!(my_subs, neighbor_subs)
    neighbor_subs[i] = left_idx
    left_rank = subs2ind(dims, neighbor_subs)

    neighbor_subs[i] = right_idx
    right_rank = subs2ind(dims, neighbor_subs)

    peer_nums[1, i] = left_rank
    peer_nums[2, i] = right_rank
  end

  return peer_nums
end

"""
  This function calculates the number of processes along each dimension
  of the grid using a factorization algorithm.  The algorithm is
  agnostic to the number of dimensions.

  I believe the algorithm produces the optimal decomposition of the processes
  into an N dimesnional grid, although I cannot formally prove it does.
  Comparison to n exhaustive search of the 1 to 500 processes in 3 dimensions
  shows identical results.

"""
function getCartesianDecomposition(comm_size::Integer, N::Integer)

  if comm_size != 1
    factor_dict = factor(comm_size)
  else
    factor_dict = Dict{Int, Int}(1 => 1)
  end
#  println("factor_dict = ", factor_dict)
  factors = collect(keys(factor_dict))
  multiplicity = collect(values(factor_dict))

  nfactors = sum(multiplicity)
#  println("factors = ", factors)
#  println("multiplicty = ", multiplicity)
#  println("nfactors = ", nfactors)

  # put factors into array, including multiples
  factor_arr = Array(Int, nfactors)
  idx = 1
  idx_cnt = 0
  for i=1:nfactors
    factor_arr[i] = factors[idx]
    idx_cnt += 1
    if idx_cnt == multiplicity[idx]
      idx += 1
      idx_cnt = 0
    end
  end

  # put N largest factors in decomp (filling with ones if not enough factors)
  sort!(factor_arr, rev=true)
  cart_decomp = Array(Int, N)
  fill!(cart_decomp, 1)
  
#  println("factor_arr = ", factor_arr)

  for i=1:min(nfactors, N)
    cart_decomp[i] = factor_arr[i]
  end

#  println("after initial insertions cart_decomp = ", cart_decomp)

  if nfactors > N
    opt_val = comm_size^(1/N)
#    println("more factors than dimensions")
    for i=(N+1):nfactors
      factor_i = factor_arr[i]
      # figure out where to put this factor in order to minimize error
      idx_opt = getOptimalIdx(cart_decomp, factor_i, opt_val)
      cart_decomp[idx_opt] *= factor_i
#      println("multiplying index ", idx_opt, " by factor ", factor_i)
#      println("cart_decomp = ", cart_decomp)

    end
  end

#  println("cart_decomp = ", cart_decomp)

  @assert prod(cart_decomp) == comm_size

  return cart_decomp
end

function getOptimalIdx(vals::AbstractArray, factor, opt_val)
# see where to add factor to vals to produce the minimum error

  min_err = typemax(Float64)
  min_idx = 0
  for i=1:length(vals)
    vals[i] *= factor
    err_i = getLSError(vals, opt_val)
    if err_i < min_err
      min_err = err_i
      min_idx = i
    end
    vals[i] = div(vals[i], factor)
  end

  return min_idx
end

function getLSError(vals::AbstractArray, opt_val)

  err = 0.0
  for i=1:length(vals)
    err_i = vals[i] - opt_val
    err += err_i*err_i
  end

  return err
end





function resize_arr(arr::AbstractMatrix)

  m, n = size(arr)
  arr2 = Array(Int, 2*m, n)

  for i=1:m
    for j=1:n
      arr2[i, j] = arr[i, j]
    end
  end

  return arr2
end








