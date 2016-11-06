include("generate_input.jl")

arr = ["10", "1.0", "1", "1"]

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

if comm_rank == 0
  if isfile("convergence.dat")
    rm("convergence.dat")
  end
end

N = 10
for i=1:5
  arr[1] = string(N)
  MPI.Barrier(MPI.COMM_WORLD)
  if comm_rank == 0
    makeinput(arr)
  end
  MPI.Barrier(MPI.COMM_WORLD)
  runcase("input.txt")
  N *= 2
end

MPI.Finalize()
if comm_rank == 0
  include("calcratio.jl")
end
