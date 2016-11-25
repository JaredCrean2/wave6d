include("generate_input.jl")
# npoints, tmax, maxdim, writeconv, Nblock, blocksize
arr = ["10", "1.0", "6", "0", "3", "2"]

#npoints = [4, 4, 4, 4, 4, 4] + 12
#tmax = [1.872535253073763e-5,1.872535253073763e-5,1.872535253073763e-5,1.872535253073763e-5,1.872535253073763e-5,1.872535253073763e-5]./10000
#tmax = zeros(length(npoints))

npoints = [16777216*4, 4096, 256, 64, 32, 16]
tmax = [1.872535253073763e-5, 0.03835888465921603, 0.410665706351607, 1.2466637514245213, 2.0268339700579308, 3.4906585039886586]./4

dir = @__FILE__
dir_rev = reverse(dir)
slash_idx = findfirst(dir_rev, '/')
path = dir_rev[ (slash_idx+1):end]
path = reverse(path)
push!(LOAD_PATH, path)

using WaveSolver
using MPI
MPI.Init()
comm_rank = MPI.Comm_rank(MPI.COMM_WORLD)
comm_size = MPI.Comm_size(MPI.COMM_WORLD)

if isdir("timing_data")
  println("deleting existing data")
  rm("timing_data", recursive=true)
end

mkdir("timing_data")
cd("./timing_data")

maxdim = 6
block_sizes = [2, 4, 8, 16]

# run blocked cases
for d=1:maxdim
  arr[1] = string(npoints[d])
  arr[2] = string(tmax[d])
  arr[3] = string(d)
  for nblock=1:d
    arr[5] = string(nblock)
    for blocksize in block_sizes
      arr[6] = string(blocksize)

      # dir name is dimension, number of blocked loops, block size
      println("running maxdim = ", d, ", nblock = ", nblock, ", blocksize = ", blocksize)

      dirname = string(d, "_", nblock, "_", blocksize)
      mkdir(dirname)
      cd(dirname)

      makeinput(arr)
      runcase("input.txt")
      gc()
      cd("..")
    end
  end
end

# run non blocked cases
nblock = 0
blocksize = 0

arr[5] = string(nblock)
arr[6] = string(blocksize)

for d=1:maxdim
  arr[1] = string(npoints[d])
  arr[2] = string(tmax[d])
  arr[3] = string(d)

  println("running maxdim = ", d, ", nblock = ", nblock, ", blocksize = ", blocksize)

  dirname = string(d, "_", nblock, "_", blocksize)
  mkdir(dirname)
  cd(dirname)

  makeinput(arr)
  runcase("input.txt")
  gc()
  cd("..")
end

# run hilbert cases
block_sizes_complete = zeros(Int, length(block_sizes) + 1)
block_sizes_complete[2:end] = block_sizes

nblock = -1
#blocksize = 0

arr[5] = string(nblock)
#arr[6] = string(blocksize)

for d=2:maxdim
  arr[1] = string(npoints[d])
  arr[2] = string(tmax[d])
  arr[3] = string(d)

  for blocksize in block_sizes_complete
    arr[6] = string(blocksize)


    println("running hilbert maxdim = ", d, ", blocksize = ", blocksize)

    dirname = string("h_", d, "_", blocksize)
    mkdir(dirname)
    cd(dirname)

    makeinput(arr)
    runcase("input.txt")
    gc()
    cd("..")
  end
end


cd("..")


