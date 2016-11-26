maxdim = 6
block_sizes = [2, 4, 8, 16]
blocksizes_extended = vcat([0], block_sizes)

nblocksizes = length(block_sizes)
nblocksizes_extended = length(blocksizes_extended)

offset = 1
# each inner array is dim + 1 by length(block_sizes)
# the 1 extra data point is  unblocked case, and the next
# dim are the number of blocked loops
data = Array(Array{Float64, 2}, maxdim)
# hilbert data
# outer array is for dimension
# inner array is for block size
hdata = Array(Array{Float64, 1}, maxdim)
for i=1:maxdim
  data[i] = Array(Float64, i + offset, nblocksizes)
  hdata[i] = Array(Float64, nblocksizes_extended)
  fill!(data[i], -1)  # use this to filter out unused entries later
  fill!(hdata[i], -1)
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
    data[d][1, blocksize_idx] = data_i[2]
  end
end

nblock = -1
blocksize = 0

for d=2:maxdim
  for blocksize_idx = 1:nblocksizes_extended
    blocksize = blocksizes_extended[blocksize_idx]

    println("getting hilbert maxdim = ", d, ", blocksize = ", blocksize)
    dirname = string("h_", d, "_", blocksize)
    fname_i = joinpath(dirname, "timing.dat")
    data_i = readdlm(fname_i)

    hdata[d][blocksize_idx] = data_i[2]

#    # use hilbert case as data point on all block sizes
#    for blocksize_idx = 1:nblocksizes
#      data[d][1, blocksize_idx] = data_i[2]
#    end
  end
end

# get data for bar chart: 0 block time, n block time, hilbert time
bardata = Array(Float64, 3, maxdim)
fill!(bardata, -1)
# nblocked loops, blocksize for the minimum time blocked loop
min_nblock_data = Array(Int, 2, maxdim)
fill!(min_nblock_data, -1)

for d=1:maxdim
  data_extract = data[d][(offset+1):end, :]
  bardata[1, d] = data[d][1, 1]
  bardata[2, d] = minimum(data_extract)

  if d == 1
    bardata[3, d] = bardata[2, d]
  else
    bardata[3, d] = minimum(hdata[d])
  end
  min_idx = findfirst(data_extract, bardata[2, d])
  min_block, min_blocksize_idx = ind2sub(data_extract, min_idx)
  min_blocksize = block_sizes[min_blocksize_idx]
  min_nblock_data[1, d] = min_block  # number of blocked loops
  min_nblock_data[2, d] = min_blocksize  # blocksize
end

# write output
for d=1:maxdim
  fname_i = string("d", d, "data.dat")
  writedlm(fname_i, data[d])
end

# write hilbert data
for d=2:maxdim
  fname_i = string("hd", d, "data.dat")
  writedlm(fname_i, hdata[d])
end

for d=1:maxdim
  fname_i = string("bar", d, "data.dat")
  writedlm(fname_i, bardata[:, d])

  fname_i = string("bar", d, "data2.dat")
  writedlm(fname_i, min_nblock_data[:, d])
end




