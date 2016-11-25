maxdim = 6
block_sizes = [2, 4, 8, 16]

nblocksizes = length(block_sizes)
offset = 2
# each inner array is dim + 2 by length(block_sizes)
# the 2 data points are the hilbert and unblocked cases, and the next
# dim are the number of blocked loops
data = Array(Array{Float64, 2}, maxdim)
for i=1:maxdim
  data[i] = Array(Float64, i + offset, nblocksizes)
  fill!(data[i], -1)  # use this to filter out unused entries later
end

for d=1:maxdim
  for nblock=1:d
    for blocksize_idx = 1:nblocksizes
       blocksize = block_sizes[blocksize_idx]

      # dir name is dimension, number of blocked loops, block size
      println("getting maxdim = ", d, ", nblock = ", nblock, ", blocksize = ", blocksize)

      dirname = string(d, "_", nblock, "_", blocksize)
      fname_i = joinpath(dirname, "timing.dat")
      data_i = readdlm(fname_i)
      # get the compute timing
      data[d][nblock + offset, blocksize_idx] = data_i[2]
    end
  end
end


# now get unblocked cases
nblock = 0
blocksize = 0

for d=1:maxdim

  println("getting unblocked dimension ", d)
  dirname = string(d, "_", nblock, "_", blocksize)
  fname_i = joinpath(dirname, "timing.dat")
  data_i = readdlm(fname_i)
  # get the compute timing
  # use the unblocked case as a data point on all block sizes
  for blocksize_idx = 1:nblocksizes
    data[d][2, blocksize_idx] = data_i[2]
  end
end
#=
nblock = -1
blocksize = 0

for d=2:maxdim

  println("getting hilbert dimension ", d)
  dirname = string("h_", d)
  fname_i = joinpath(dirname, "timing.dat")
  data_i = readdlm(fname_i)
  # use hilbert case as data poitn on all block sizes
  for blocksize_idx = 1:nblocksizes
    data[d][1, blocksize_idx] = data_i[2]
  end
end

# a 1d hilbert curve is the unblocked case
data[1][1, :] = data[1][2, :]  
=#
# write output
for d=1:maxdim
  fname_i = string("d", d, "data.dat")
  writedlm(fname_i, data[d])
end




