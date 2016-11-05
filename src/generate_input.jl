
if length(ARGS) != 3
  println("usage: N, tmax, maxdim")
  exit()
end

N = ARGS[1]
tmax = ARGS[2]
maxdim = parse(Int, ARGS[3])


f = open("input.txt", "w")
for i=1:maxdim
  println(f, N)
end


for i=1:maxdim
  println(f, 0.0)
  println(f, 2*pi)
end

println(f, tmax)

close(f)
