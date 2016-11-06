# function to generate input from array of strings

function makeinput(arr)

  if length(arr) != 4
    println("usage: N, tmax, maxdim, convergence")
    exit()
  end

  N = arr[1]
  tmax = arr[2]
  maxdim = parse(Int, arr[3])
  convergence = parse(Int, arr[4])

  f = open("input.txt", "w")
  for i=1:maxdim
    println(f, N)
  end


  for i=1:maxdim
    println(f, 0.0)
    println(f, 2*pi)
  end

  println(f, tmax)
  println(f, convergence)

  close(f)
end
