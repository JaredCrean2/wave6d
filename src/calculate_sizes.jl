# calculate the number of points and max time for each dimension given
# the size of the array and number of time steps

xmin = 0.0
xmax = 2*pi
arr_size = 128.0 # size in megabytes
nsteps = 100

function getDeltaT(delta_xs, CFL)

  val = 0.0
  for i=1:length(delta_xs)
    val += 1/delta_xs[i]
  end

  return CFL/val
end

function getNearestMultiple(n, fac)
# get the nearest multiple of fac to n

  if n % fac == 0
    return n
  end

  if (n/fac - div(n, fac)) > 0.5
    for i=1:fac
      n += 1
      if (n % fac) == 0
        return n
      end
    end
  else
    for i=1:fac
      n -= 1
      if (n % fac) == 0
        return n
      end
    end
  end

  return n
end


for d=1:6
  delta_xs = zeros(d)
  arr_size_i = 0.0
  n = 1
  while (arr_size_i < arr_size)
    n += 1
    real_n = n + 4  # include ghost points
    arr_size_i = ((real_n^d)*8)/(1024*1024)
  end

  # calculate delta_xs
  fill!(delta_xs, (xmax - xmin)/(n-1))
  delta_t = getDeltaT(delta_xs, 0.5)
#  println("delta_t = ", delta_t)

  tmax = 100*delta_t
  n_nearest = getNearestMultiple(n, 8)
  println("for dimension ", d, ", number of points = ", n, ", array size = ", arr_size_i, ", tmax = ", tmax)
  println("nearest multiple = ", n_nearest)
end




