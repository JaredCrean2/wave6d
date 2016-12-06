# run all convergence studies

include("generate_input.jl")
# npoints, tmax, maxdim, writeconv, Nblock, blocksize
arr = ["10", "1.0", "6", "1", "3", "2"]

function runconv(N, other_vals)

  for i=1:4
    arr[1] = string(N)
    MPI.Barrier(MPI.COMM_WORLD)
    if comm_rank == 0
      makeinput(arr)
    end
    MPI.Barrier(MPI.COMM_WORLD)
    runcase("input.txt")
    N *= 2
  end
end




#npoints = [4, 4, 4, 4, 4, -12] + 28
#tmax = [1.872535253073763e-5,1.872535253073763e-5,1.872535253073763e-5,1.872535253073763e-5,1.872535253073763e-5,1.872535253073763e-5]./10000
#tmax = zeros(length(npoints))


npoints = [4, 4, 4, 4, 4, 4 - 16] + 28
tmax = [1.0, 0.5, 0.5, 0.25, 0.25, 0.25]
#tmax = zeros(length(npoints))
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


if isdir("convergence_data")
  println("deleting existing data")
  rm("convergence_data", recursive=true)
end


mkdir("convergence_data")

cd("./convergence_data")

maxdim = 6
block_sizes = [4, 8, 16, 32]

# run blocked cases
for d=1:maxdim
  arr[1] = string(npoints[d])
  arr[2] = string(tmax[d])
  arr[3] = string(d)
  for nblock=1:d
    arr[5] = string(nblock)
    for blocksize in block_sizes
      arr[6] = string(blocksize)

      # skip really large case
      if d == 6 && blocksize > 16
        continue
      end

      # dir name is dimension, number of blocked loops, block size
      println("running maxdim = ", d, ", nblock = ", nblock, ", blocksize = ", blocksize)

      dirname = string(d, "_", nblock, "_", blocksize)
      mkdir(dirname)
      cd(dirname)

#      makeinput(arr)
#      runcase("input.txt")
       runconv(npoints[d], arr)
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

#  makeinput(arr)
#  runcase("input.txt")

  runconv(npoints[d], arr)
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

    # skip really large case
    if d == 6 && blocksize > 16
      continue
    end

    println("running hilbert maxdim = ", d, ", blocksize = ", blocksize)

    dirname = string("h_", d, "_", blocksize)
    mkdir(dirname)
    cd(dirname)

#    makeinput(arr)
#    runcase("input.txt")   
    runconv(npoints[d], arr)
    gc()
    cd("..")
  end
end


cd("..")


