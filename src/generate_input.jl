# function to generate input from array of strings

function makeinput(arr)

  if length(arr) < 4 || length(arr) > 6
    println("usage: N, tmax, maxdim, convergence Nblock=0 blocksize=4")
    exit()
  end

  N = arr[1]
  tmax = arr[2]
  maxdim = parse(Int, arr[3])
  convergence = parse(Int, arr[4])
  if length(arr) >= 5
    Nblock = parse(Int, arr[5])
  else
    Nblock = 0
  end

  println("Nblock = ", Nblock)

  if length(arr) >= 6
    blocksize = parse(Int, arr[6])
  else
    blocksize = 4
  end

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

  println(f, Nblock)
  println(f, blocksize)

  close(f)
end
