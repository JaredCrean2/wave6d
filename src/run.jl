if length(ARGS) != 1
  println("usage: julia run.jl name_of_input_file")
  exit(1)
end

dir = @__FILE__
dir_rev = reverse(dir)
slash_idx = findfirst(dir_rev, '/')
path = dir_rev[ (slash_idx+1):end]
path = reverse(path)
push!(LOAD_PATH, path)

using WaveSolver



MPI.Init()
runcase(ARGS[1])

MPI.Finalize()
