# interface functions for hilbert curve needed by the solver

typealias bitmask_t Culonglong

# get the arrays needed for a given dimension N
# npoints is the number of points on the curve to get at one time
function getState(N::Integer, npoints::Integer)
  coords = zeros(bitmask_t, N)
  idxs = zeros(UInt16, N, npoints)  # use UInt16s to save space
  
  return coords, idxs
end

# void loadNpoints(int n, int dims, bitmask_t coords[], uint16_t idxs[]);
function loadNpoints(npoints::Integer, ndims::Integer, 
                     coords::AbstractArray{bitmask_t, 1}, 
                     idxs::AbstractArray{UInt16, 2})

  ccall( (:loadNpoints, "hilbert"), Void, (Cint, Cint, Ptr{bitmask_t}, Ptr{UInt16}), npoints, ndims, coords, idxs)

  return nothing
end


#int checkDimensions(int dims, int npoints);
function checkDimensions(maxdim::Integer, max_points::Integer)

  ret_status = ccall( (:checkDimensions, "hilbert"), Cint, (Cint, Cint), maxdim, max_points);

  if ret_status != 0
    throw(ErrorException("checkDimensions returned non-zero exit status"))
  end

  return nothing
end

function getNumBlocks(num_grid_points::Integer, maxdim::Integer, blocksize::Integer)
# for a grid with num_grid_points in each direction and maxdim directions,
# calculate how many blocks of blocksize points it takes to traverse the
# entire grid

  num_grid_points_big = Int128(num_grid_points)
  total_grid_points = num_grid_points_big

  for i=2:maxdim
    total_grid_points *= num_grid_points_big
  end

  if total_grid_points % blocksize != 0
    throw(ErrorException("grid of $num_grid_points_big is not divisible by a blocksize of $blocksize"))
  end

  return div(total_grid_points, blocksize)
end

  
