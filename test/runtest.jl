using FactCheck

push!(LOAD_PATH, string(pwd(), "/../src"))
using WaveSolver

function hasrow(arr, row)

  for i=1:size(arr, 1)
    if arr[i, :] == row
      return true
    end
  end

  return false
end



# test the calculation of the MPI rank grid
#=
facts("----- Testing getMPIMatches -----") do
  matches, idx = WaveSolver.getMPIMatches(4, 2)
  println("typeof(matches) = ", typeof(matches))

  println("matches = \n", matches[1:idx, :])


  @fact size(matches, 2) --> 2
  @fact hasrow(matches, [1 4]) --> true
  @fact hasrow(matches, [4 1]) --> true
  @fact hasrow(matches, [2 2]) --> true

  matches, idx = WaveSolver.getMPIMatches(27, 3)
  println("mataches = \n", matches[1:idx, :])
  @fact hasrow(matches, [3  3 3]) --> true

end
=#

function testBig(comm_size)

  N = 6
  WaveSolver.getMPIMatches(comm_size, N)
end

function exhaustiveSearch(comm_size, N)

  matches, idx = WaveSolver.getMPIMatches(comm_size, N)

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


  return vec(matches[idx_opt, :])
end

function checkPermutiveEquality(x1, x2)

  @assert length(x1) == length(x2)

  allmatch = true
  for i=1:length(x1)
    val = x1[i]
    pred = x -> x == val
    c1 = count(pred, x1)
    c2 = count(pred, x2)

    allmatch = allmatch && (c1 == c2)
  end

  return allmatch
end
 

#=
@time testBig(10)
@time testBig(20)
@time testBig(30)
@time testBig(40)
@time testBig(100)

=#

facts("----- testing getCartesianDecomposition -----") do
  cart = WaveSolver.getCartesianDecomposition(1, 1)
  println("cart = ", cart, "\n")
  @fact cart --> [1]
  cart = WaveSolver.getCartesianDecomposition(2, 1)
  @fact cart --> [2]
  println("\ncart 2, 1 = ", cart)

  cart = WaveSolver.getCartesianDecomposition(2, 2)
  @fact cart --> [2, 1]
  println("cart 2 2 = ", cart)

  cart = WaveSolver.getCartesianDecomposition(3, 2)
  @fact cart --> [3, 1]
 
  cart = WaveSolver.getCartesianDecomposition(4, 2)
  @fact cart --> [2, 2]

  cart = WaveSolver.getCartesianDecomposition(5, 2)
  @fact cart --> [5, 1]

  # 3d
  cart = WaveSolver.getCartesianDecomposition(2, 3)
  @fact cart --> [2, 1, 1]

  print("\n")
  cart = WaveSolver.getCartesianDecomposition(2*3*5, 3)
  @fact cart --> [5, 3, 2]
 
  print("\n")
  cart = WaveSolver.getCartesianDecomposition(2*3*4*5, 3)
  @fact cart --> [5, 6, 4]
 
  cart = WaveSolver.getCartesianDecomposition(2*3*5*7*11*13, 3)
  @fact cart --> [26, 33, 35]

  x1 = [1 2 3]
  x2 = [2 3 1]
  
  @fact checkPermutiveEquality(x1, x2) --> true

  # check against exhaustive search
  @time for i=1:64
#    println("i = ", i)
    cart1 = WaveSolver.getCartesianDecomposition(i, 3)
    cart2 = exhaustiveSearch(i, 3)
    @fact checkPermutiveEquality(x1, x2) --> true
  end
 
end

facts("----- testing getGridInfo -----") do


  dim_vec = WaveSolver.getCartesianDecomposition(1, 1)
  peer_nums, my_subs = WaveSolver.getGridInfo(dim_vec, 0)
  peer_nums += 1
  @fact peer_nums --> ones(2, 1)
  @fact my_subs --> [1]

  dim_vec = WaveSolver.getCartesianDecomposition(2, 1)
  peer_nums, my_subs = WaveSolver.getGridInfo(dim_vec, 0)
  peer_nums += 1

  exp_peernums = 2*ones(2, 1)
  @fact peer_nums --> exp_peernums
  @fact my_subs --> [1]

  peer_nums, my_subs = WaveSolver.getGridInfo(dim_vec, 1)
  peer_nums += 1
  @fact peer_nums --> ones(2, 1)
  @fact my_subs --> [2]

  println("\ntesting 2d GridInfo")
  myrank = 0
  dim_vec = WaveSolver.getCartesianDecomposition(4, 2)
  peer_nums, my_subs = WaveSolver.getGridInfo(dim_vec, myrank)
  peer_nums += 1
  println("dim_vec = ", dim_vec)
  println("peer_nums = ", peer_nums)
  println("my_subs = ", my_subs)
  exp_peer_nums = zeros(Int, 2, 2)
  exp_peer_nums[:, 1] = [2, 2]
  exp_peer_nums[:, 2] = [3, 3]
  @fact peer_nums --> exp_peer_nums
  @fact my_subs --> [1, 1]

  myrank = 1
  peer_nums, my_subs = WaveSolver.getGridInfo(dim_vec, myrank)
  peer_nums += 1
  exp_peer_nums[:, 1] = [1, 1]
  exp_peer_nums[:, 2] = [4, 4]
  @fact peer_nums --> exp_peer_nums
  @fact my_subs --> [2, 1]

  println("\ntesting myrank = 2")
  myrank = 2
  peer_nums, my_subs = WaveSolver.getGridInfo(dim_vec, myrank)
  peer_nums += 1
  exp_peer_nums[:, 1] = [4, 4]
  exp_peer_nums[:, 2] = [1, 1]
  @fact peer_nums --> exp_peer_nums
  @fact my_subs --> [1, 2]

  myrank = 3
  peer_nums, my_subs = WaveSolver.getGridInfo(dim_vec, myrank)
  peer_nums += 1
  exp_peer_nums[:, 1] = [3, 3]
  exp_peer_nums[:, 2] = [2, 2]
  @fact peer_nums --> exp_peer_nums
  @fact my_subs --> [2, 2]

end

facts("----- testing getNumPoints -----") do

  nghost = 2
  my_subs = [1]
  dim_vec = [1]
  npoints = [20]
  local_points, global_points = WaveSolver.getNumPoints(my_subs, dim_vec, npoints, nghost)

  npoints_global_exp = zeros(Int, 2, 1)
  npoints_global_exp[1] = 1
  npoints_global_exp[2] = 20
  @fact local_points --> npoints
  @fact global_points --> npoints_global_exp

  # test with 2 processes
  my_subs = [1]
  dim_vec[1] = 2
  npoints = [20]

  local_points, global_points = WaveSolver.getNumPoints(my_subs, dim_vec, npoints, nghost)
  npoints_global_exp[1] = 1
  npoints_global_exp[2] = 10

  @fact local_points --> [10]
  @fact global_points --> npoints_global_exp

  my_subs[1] = 2

  local_points, global_points = WaveSolver.getNumPoints(my_subs, dim_vec, npoints, nghost)
  npoints_global_exp[1] = 11
  npoints_global_exp[2] = 20

  @fact local_points --> [10]
  @fact global_points --> npoints_global_exp


  # test with 4 processes in 4 dimensions
  my_subs = [1, 1]
  dim_vec = [2, 2]
  npoints = [10, 10]
  npoints_global_exp = zeros(Int, 2, 2)
  npoints_local_exp = [5, 5]

  npoints_global_exp[:, 1] = [1, 5]
  npoints_global_exp[:, 2] = [1, 5]

  local_points, global_points = WaveSolver.getNumPoints(my_subs, dim_vec, npoints, nghost)

  @fact local_points --> npoints_local_exp
  @fact global_points --> npoints_global_exp

  my_subs = [2, 1]
  npoints_global_exp[:, 1] = [6, 10]
  npoints_global_exp[:, 2] = [1, 5]

  local_points, global_points = WaveSolver.getNumPoints(my_subs, dim_vec, npoints, nghost)

  @fact local_points --> npoints_local_exp
  @fact global_points --> npoints_global_exp

  my_subs = [1, 2]
  npoints_global_exp[:, 1] = [1, 5]
  npoints_global_exp[:, 2] = [6, 10]

  local_points, global_points = WaveSolver.getNumPoints(my_subs, dim_vec, npoints, nghost)

  @fact local_points --> npoints_local_exp
  @fact global_points --> npoints_global_exp

  my_subs = [2, 2]
  npoints_global_exp[:, 1] = [6, 10]
  npoints_global_exp[:, 2] = [6, 10]

  local_points, global_points = WaveSolver.getNumPoints(my_subs, dim_vec, npoints, nghost)

  @fact local_points --> npoints_local_exp
  @fact global_points --> npoints_global_exp

  # test non evenly divisible cases
  my_subs = [1]
  dim_vec = [2]
  npoints = [9]
  npoints_global_exp = zeros(Int, 2, 1)
  npoints_local_exp = [5]

  npoints_global_exp[:, 1] = [1, 5]

  local_points, global_points = WaveSolver.getNumPoints(my_subs, dim_vec, npoints, nghost)
  @fact local_points --> npoints_local_exp
  @fact global_points --> npoints_global_exp


  my_subs = [2]
  npoints_global_exp = zeros(Int, 2, 1)
  npoints_local_exp = [4]

  npoints_global_exp[:, 1] = [6, 9]

  local_points, global_points = WaveSolver.getNumPoints(my_subs, dim_vec, npoints, nghost)
  @fact local_points --> npoints_local_exp
  @fact global_points --> npoints_global_exp



end

FactCheck.exitstatus()
