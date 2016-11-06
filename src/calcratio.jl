
data = readdlm("convergence.dat")
m, n = size(data)
data2 = zeros(m, n+1)
@assert n == 2

data2[:, 1:n] = data[:, :]

for i=2:m
  n1 = 1/data2[i, 1]
  n2 = 1/data2[i-1, 1]
  data2[i, n+1] = log(data2[i-1, n]/data2[i, n])/log(n1/n2)
end

println("data2 = \n", data2)

#writedlm("convergence2.dat", data2)
f = open("convergence2.dat", "w")
println(f, "delta_x  Max. Err  Convergence Rate")
for i=1:m
  @printf(f, "%f %e %f\n", data2[i, 1], data2[i,2], data2[i,3])
end
