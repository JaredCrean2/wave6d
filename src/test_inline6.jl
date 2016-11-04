# test whether putting the equation index as the first or the last is better

@inline function kernel(u_i::AbstractArray, u_ip1::AbstractArray, idx)

  d1 = idx[1]
  d2 = idx[2]
  d3 = idx[3]
  d4 = idx[4]
  d5 = idx[5]

  delta_12 = 0.5
  delta_22 = delta_12
  delta_32 = delta_12
  delta_42 = delta_12
  delta_52 = delta_12


  u11 = (1.0)*u_i[ 1, d1 - 2, d2, d3, d4, d5] +
        (2.0)*u_i[ 1, d1 - 1, d2, d3, d4, d5] +
        (3.0)*u_i[ 1, d1    , d2, d3, d4, d5] +
        (4.0)*u_i[ 1, d1 + 1, d2, d3, d4, d5] +
        (5.0)*u_i[ 1, d1 + 2, d2, d3, d4, d5]

  u22 = (1.0)*u_i[ 1, d1, d2 - 2, d3, d4, d5] +
        (2.0)*u_i[ 1, d1, d2 - 1, d3, d4, d5] +
        (3.0)*u_i[ 1, d1, d2    , d3, d4, d5] +
        (4.0)*u_i[ 1, d1, d2 + 1, d3, d4, d5] +
        (5.0)*u_i[ 1, d1, d2 + 2, d3, d4, d5]

  u33 = (1.0)*u_i[ 1, d1, d2, d3 - 2, d4, d5] +
        (2.0)*u_i[ 1, d1, d2, d3 - 1, d4, d5] +
        (3.0)*u_i[ 1, d1, d2, d3    , d4, d5] +
        (4.0)*u_i[ 1, d1, d2, d3 + 1, d4, d5] +
        (5.0)*u_i[ 1, d1, d2, d3 + 2, d4, d5]

  u44 = (1.0)*u_i[ 1, d1, d2, d3, d4 - 2, d5] +
        (2.0)*u_i[ 1, d1, d2, d3, d4 - 1, d5] +
        (3.0)*u_i[ 1, d1, d2, d3, d4    , d5] +
        (4.0)*u_i[ 1, d1, d2, d3, d4 + 1, d5] +
        (5.0)*u_i[ 1, d1, d2, d3, d4 + 2, d5]

  u55 = (1.0)*u_i[ 1, d1, d2, d3, d4, d5 - 2] +
        (2.0)*u_i[ 1, d1, d2, d3, d4, d5 - 1] +
        (3.0)*u_i[ 1, d1, d2, d3, d4, d5    ] +
        (4.0)*u_i[ 1, d1, d2, d3, d4, d5 + 1] +
        (5.0)*u_i[ 1, d1, d2, d3, d4, d5 + 2]

  u_ip1[2,  d1, d2, d3, d4, d5] = delta_12*u11 + delta_22*u22 + delta_32*u33 + 
                               delta_42*u44 + delta_52*u55

  u11 = (1.0)*u_i[ 2, d1 - 2, d2, d3, d4, d5] +
        (2.0)*u_i[ 2, d1 - 1, d2, d3, d4, d5] +
        (3.0)*u_i[ 2, d1    , d2, d3, d4, d5] +
        (4.0)*u_i[ 2, d1 + 1, d2, d3, d4, d5] +
        (5.0)*u_i[ 2, d1 + 2, d2, d3, d4, d5]

  u22 = (1.0)*u_i[ 2, d1, d2 - 2, d3, d4, d5] +
        (2.0)*u_i[ 2, d1, d2 - 1, d3, d4, d5] +
        (3.0)*u_i[ 2, d1, d2    , d3, d4, d5] +
        (4.0)*u_i[ 2, d1, d2 + 1, d3, d4, d5] +
        (5.0)*u_i[ 2, d1, d2 + 2, d3, d4, d5]

  u33 = (1.0)*u_i[ 2, d1, d2, d3 - 2, d4, d5] +
        (2.0)*u_i[ 2, d1, d2, d3 - 1, d4, d5] +
        (3.0)*u_i[ 2, d1, d2, d3    , d4, d5] +
        (4.0)*u_i[ 2, d1, d2, d3 + 1, d4, d5] +
        (5.0)*u_i[ 2, d1, d2, d3 + 2, d4, d5]

  u44 = (1.0)*u_i[ 2, d1, d2, d3, d4 - 2, d5] +
        (2.0)*u_i[ 2, d1, d2, d3, d4 - 1, d5] +
        (3.0)*u_i[ 2, d1, d2, d3, d4    , d5] +
        (4.0)*u_i[ 2, d1, d2, d3, d4 + 1, d5] +
        (5.0)*u_i[ 2, d1, d2, d3, d4 + 2, d5]

  u55 = (1.0)*u_i[ 2, d1, d2, d3, d4, d5 - 2] +
        (2.0)*u_i[ 2, d1, d2, d3, d4, d5 - 1] +
        (3.0)*u_i[ 2, d1, d2, d3, d4, d5    ] +
        (4.0)*u_i[ 2, d1, d2, d3, d4, d5 + 1] +
        (5.0)*u_i[ 2, d1, d2, d3, d4, d5 + 2]

  u_ip1[2, d1, d2, d3, d4, d5] = delta_12*u11 + delta_22*u22 + delta_32*u33 + 
                               delta_42*u44 + delta_52*u55






  return nothing
end


function outer_func(u_i, u_ip1)

  s1 = size(u_i, 3)
  nghost = 2
  ia = nghost + 1
  ib = s1 - nghost

    for d5=ia:ib
      for d4=ia:ib
        for d3 = ia:ib
          for d2= ia:ib
            @simd for d1 = ia:ib
              idx = (d1, d2, d3, d4, d5)
              kernel(u_i, u_ip1, idx)
            end
          end
        end
      end
    end

  return nothing
end

@inline function kernel2(u_i::AbstractArray, u_ip1::AbstractArray, idx)

  d1 = idx[1]
  d2 = idx[2]
  d3 = idx[3]
  d4 = idx[4]
  d5 = idx[5]

  delta_12 = 0.5
  delta_22 = delta_12
  delta_32 = delta_12
  delta_42 = delta_12
  delta_52 = delta_12


  u11 = (1.0)*u_i[ d1 - 2, d2, d3, d4, d5, 1] +
        (2.0)*u_i[ d1 - 1, d2, d3, d4, d5, 1] +
        (3.0)*u_i[ d1    , d2, d3, d4, d5, 1] +
        (4.0)*u_i[ d1 + 1, d2, d3, d4, d5, 1] +
        (5.0)*u_i[ d1 + 2, d2, d3, d4, d5, 1]

  u22 = (1.0)*u_i[ d1, d2 - 2, d3, d4, d5, 1] +
        (2.0)*u_i[ d1, d2 - 1, d3, d4, d5, 1] +
        (3.0)*u_i[ d1, d2    , d3, d4, d5, 1] +
        (4.0)*u_i[ d1, d2 + 1, d3, d4, d5, 1] +
        (5.0)*u_i[ d1, d2 + 2, d3, d4, d5, 1]

  u33 = (1.0)*u_i[ d1, d2, d3 - 2, d4, d5, 1] +
        (2.0)*u_i[ d1, d2, d3 - 1, d4, d5, 1] +
        (3.0)*u_i[ d1, d2, d3    , d4, d5, 1] +
        (4.0)*u_i[ d1, d2, d3 + 1, d4, d5, 1] +
        (5.0)*u_i[ d1, d2, d3 + 2, d4, d5, 1]

  u44 = (1.0)*u_i[ d1, d2, d3, d4 - 2, d5, 1] +
        (2.0)*u_i[ d1, d2, d3, d4 - 1, d5, 1] +
        (3.0)*u_i[ d1, d2, d3, d4    , d5, 1] +
        (4.0)*u_i[ d1, d2, d3, d4 + 1, d5, 1] +
        (5.0)*u_i[ d1, d2, d3, d4 + 2, d5, 1]

  u55 = (1.0)*u_i[ d1, d2, d3, d4, d5 - 2, 1] +
        (2.0)*u_i[ d1, d2, d3, d4, d5 - 1, 1] +
        (3.0)*u_i[ d1, d2, d3, d4, d5    , 1] +
        (4.0)*u_i[ d1, d2, d3, d4, d5 + 1, 1] +
        (5.0)*u_i[ d1, d2, d3, d4, d5 + 2, 1]

  u_ip1[ d1, d2, d3, d4, d5, 2] = delta_12*u11 + delta_22*u22 + delta_32*u33 + 
                               delta_42*u44 + delta_52*u55


  u11 = (1.0)*u_i[ d1 - 2, d2, d3, d4, d5, 2] +
        (2.0)*u_i[ d1 - 1, d2, d3, d4, d5, 2] +
        (3.0)*u_i[ d1    , d2, d3, d4, d5, 2] +
        (4.0)*u_i[ d1 + 1, d2, d3, d4, d5, 2] +
        (5.0)*u_i[ d1 + 2, d2, d3, d4, d5, 2]

  u22 = (1.0)*u_i[ d1, d2 - 2, d3, d4, d5, 2] +
        (2.0)*u_i[ d1, d2 - 1, d3, d4, d5, 2] +
        (3.0)*u_i[ d1, d2    , d3, d4, d5, 2] +
        (4.0)*u_i[ d1, d2 + 1, d3, d4, d5, 2] +
        (5.0)*u_i[ d1, d2 + 2, d3, d4, d5, 2]

  u33 = (1.0)*u_i[ d1, d2, d3 - 2, d4, d5, 2] +
        (2.0)*u_i[ d1, d2, d3 - 1, d4, d5, 2] +
        (3.0)*u_i[ d1, d2, d3    , d4, d5, 2] +
        (4.0)*u_i[ d1, d2, d3 + 1, d4, d5, 2] +
        (5.0)*u_i[ d1, d2, d3 + 2, d4, d5, 2]

  u44 = (1.0)*u_i[ d1, d2, d3, d4 - 2, d5, 2] +
        (2.0)*u_i[ d1, d2, d3, d4 - 1, d5, 2] +
        (3.0)*u_i[ d1, d2, d3, d4    , d5, 2] +
        (4.0)*u_i[ d1, d2, d3, d4 + 1, d5, 2] +
        (5.0)*u_i[ d1, d2, d3, d4 + 2, d5, 2]

  u55 = (1.0)*u_i[ d1, d2, d3, d4, d5 - 2, 2] +
        (2.0)*u_i[ d1, d2, d3, d4, d5 - 1, 2] +
        (3.0)*u_i[ d1, d2, d3, d4, d5    , 2] +
        (4.0)*u_i[ d1, d2, d3, d4, d5 + 1, 2] +
        (5.0)*u_i[ d1, d2, d3, d4, d5 + 2, 2]

  u_ip1[ d1, d2, d3, d4, d5, 1] = delta_12*u11 + delta_22*u22 + delta_32*u33 + 
                               delta_42*u44 + delta_52*u55



  return nothing
end


function outer_func2(u_i, u_ip1)

  s1 = size(u_i, 3)
  nghost = 2
  ia = nghost + 1
  ib = s1 - nghost

    for d5=ia:ib
      for d4=ia:ib
        for d3 = ia:ib
          for d2= ia:ib
            @simd for d1 = ia:ib
              idx = (d1, d2, d3, d4, d5)
              kernel2(u_i, u_ip1, idx)
            end
          end
        end
      end
    end

  return nothing
end



function runtest()
  # dim = 20 is used for benchmarking
  dim = 30
  arr_size = (dim^6)*8/(1024^2)
  println("array size = ", arr_size, " megabytes")

  arr1 = rand(2, dim, dim, dim, dim, dim)
  arr2 = zeros(size(arr1))


  println("testing outer_func")
  @time outer_func(arr1, arr2)
  fill!(arr2, 0.0)
  @time outer_func(arr1, arr2)
  fill!(arr2, 0.0)



  arr1 = rand( dim, dim, dim, dim, dim, 2)
  arr2 = zeros(size(arr1))


  println("testing outer_func2")

  @time outer_func2(arr1, arr2)
  fill!(arr2, 0.0)
  @time outer_func2(arr1, arr2)
  fill!(arr2, 0.0)

#=
  dtype = typeof(arr1)

  fname = "outer_func_llvm.txt"
  f = open(fname, "w")
  println(f, "outer_func code_llvm = ")
  code_llvm(f, outer_func, (dtype, dtype))
  close(f)

  fname = "outer_func_native.txt"
  f = open(fname, "w")
  println(f, "outer_func code_native = ")
  code_native(f, outer_func, (dtype, dtype))
  close(f)




  fname = "outer_func2_llvm.txt"
  f = open(fname, "w")
  println(f, "outer_func2 code_llvm = ")
  code_llvm(f, outer_func2, (dtype, dtype))
  close(f)

  fname = "outer_func2_native.txt"
  f = open(fname, "w")
  println(f, "outer_func2 code_native = ")
  code_native(f, outer_func2, (dtype, dtype))
  close(f)

=#

end

runtest()


