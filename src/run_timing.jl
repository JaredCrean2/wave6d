include("generate_input.jl")
# npoints, tmax, maxdim, writeconv, Nblock, blocksize
arr = ["10", "1.0", "6", "0", "3", "2"]

#npoints = [4, 4, 4, 4, 4, 4] + 4
npoints = [16777208, 4088, 248, 56, 24, 16]
tmax = [1.8725356995210904e-5, 0.03839639029075767, 0.4172101797596007, 1.331183327792285, 2.731819698773733, 4.759988869075444]

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
  rm("timing_data", recursive=true)
end
mkdir("timing_data")
cd("./timing_data")

maxdim = 6
block_sizes = [2, 4, 8]

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
  cd("..")
end

cd("..")


