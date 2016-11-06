include("generate_input.jl")

arr = ["10", "1.0", "1", "1"]

dir = @__FILE__
dir_rev = reverse(dir)
slash_idx = findfirst(dir_rev, '/')
path = dir_rev[ (slash_idx+1):end]
path = reverse(path)
push!(LOAD_PATH, path)

using WaveSolver
MPI.Init()

N = 10
for i=1:5
  arr[1] = string(N)
  makeinput(arr)
  runcase("input.txt")
  N *= 2
end

MPI.Finalize()
