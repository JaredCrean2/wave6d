
@inline function kernel(u_i::AbstractArray, u_ip1::AbstractArray, idx)

  d1 = idx[1]
  d2 = idx[2]
  d3 = idx[3]
  d4 = idx[4]
  d5 = idx[5]
  d6 = idx[6]

  delta_12 = 0.5
  delta_22 = delta_12
  delta_32 = delta_12
  delta_42 = delta_12
  delta_52 = delta_12
  delta_62 = delta_12
 
  u11 = (1.0)*u_i[ d1 - 2, d2, d3, d4, d5, d6] +
        (2.0)*u_i[ d1 - 1, d2, d3, d4, d5, d6] +
        (3.0)*u_i[ d1    , d2, d3, d4, d5, d6] +
        (4.0)*u_i[ d1 + 1, d2, d3, d4, d5, d6] +
        (5.0)*u_i[ d1 + 2, d2, d3, d4, d5, d6]

  u22 = (1.0)*u_i[ d1, d2 - 2, d3, d4, d5, d6] +
        (2.0)*u_i[ d1, d2 - 1, d3, d4, d5, d6] +
        (3.0)*u_i[ d1, d2    , d3, d4, d5, d6] +
        (4.0)*u_i[ d1, d2 + 1, d3, d4, d5, d6] +
        (5.0)*u_i[ d1, d2 + 2, d3, d4, d5, d6]

  u33 = (1.0)*u_i[ d1, d2, d3 - 2, d4, d5, d6] +
        (2.0)*u_i[ d1, d2, d3 - 1, d4, d5, d6] +
        (3.0)*u_i[ d1, d2, d3    , d4, d5, d6] +
        (4.0)*u_i[ d1, d2, d3 + 1, d4, d5, d6] +
        (5.0)*u_i[ d1, d2, d3 + 2, d4, d5, d6]

  u44 = (1.0)*u_i[ d1, d2, d3, d4 - 2, d5, d6] +
        (2.0)*u_i[ d1, d2, d3, d4 - 1, d5, d6] +
        (3.0)*u_i[ d1, d2, d3, d4    , d5, d6] +
        (4.0)*u_i[ d1, d2, d3, d4 + 1, d5, d6] +
        (5.0)*u_i[ d1, d2, d3, d4 + 2, d5, d6]

  u55 = (1.0)*u_i[ d1, d2, d3, d4, d5 - 2, d6] +
        (2.0)*u_i[ d1, d2, d3, d4, d5 - 1, d6] +
        (3.0)*u_i[ d1, d2, d3, d4, d5    , d6] +
        (4.0)*u_i[ d1, d2, d3, d4, d5 + 1, d6] +
        (5.0)*u_i[ d1, d2, d3, d4, d5 + 2, d6]

  u66 = (1.0)*u_i[ d1, d2, d3, d4, d5, d6 - 2] +
        (2.0)*u_i[ d1, d2, d3, d4, d5, d6 - 1] +
        (3.0)*u_i[ d1, d2, d3, d4, d5, d6    ] +
        (4.0)*u_i[ d1, d2, d3, d4, d5, d6 + 1] +
        (5.0)*u_i[ d1, d2, d3, d4, d5, d6 + 2]

  u_ip1[ d1, d2, d3, d4, d5, d6] = delta_12*u11 + delta_22*u22 + delta_32*u33 + 
                                   delta_42*u44 + delta_52*u55 + delta_62*u66



  return nothing
end


function outer_func(u_i, u_ip1)

  s1 = size(u_i, 1)
  nghost = 2
  ia = nghost + 1
  ib = s1 - nghost

  for d1=ia:ib
    for d2=ia:ib
      for d3=ia:ib
        for d4 = ia:ib
          for d5= ia:ib
            for d6 = ia:ib
              idx = (d1, d2, d3, d4, d5, d5)
              kernel(u_i, u_ip1, idx)
            end
          end
        end
      end
    end
  end

  return nothing
end


function outer_func2(u_i, u_ip1)
  s1 = size(u_i, 1)
  nghost = 2
  ia = nghost + 1
  ib = s1 - nghost

  for d1=ia:ib
    for d2=ia:ib
      for d3=ia:ib
        for d4 = ia:ib
          for d5= ia:ib
            for d6 = ia:ib


              delta_12 = 0.5
              delta_22 = delta_12
              delta_32 = delta_12
              delta_42 = delta_12
              delta_52 = delta_12
              delta_62 = delta_12
             
              u11 = (1.0)*u_i[ d1 - 2, d2, d3, d4, d5, d6] +
                    (2.0)*u_i[ d1 - 1, d2, d3, d4, d5, d6] +
                    (3.0)*u_i[ d1    , d2, d3, d4, d5, d6] +
                    (4.0)*u_i[ d1 + 1, d2, d3, d4, d5, d6] +
                    (5.0)*u_i[ d1 + 2, d2, d3, d4, d5, d6]

              u22 = (1.0)*u_i[ d1, d2 - 2, d3, d4, d5, d6] +
                    (2.0)*u_i[ d1, d2 - 1, d3, d4, d5, d6] +
                    (3.0)*u_i[ d1, d2    , d3, d4, d5, d6] +
                    (4.0)*u_i[ d1, d2 + 1, d3, d4, d5, d6] +
                    (5.0)*u_i[ d1, d2 + 2, d3, d4, d5, d6]

              u33 = (1.0)*u_i[ d1, d2, d3 - 2, d4, d5, d6] +
                    (2.0)*u_i[ d1, d2, d3 - 1, d4, d5, d6] +
                    (3.0)*u_i[ d1, d2, d3    , d4, d5, d6] +
                    (4.0)*u_i[ d1, d2, d3 + 1, d4, d5, d6] +
                    (5.0)*u_i[ d1, d2, d3 + 2, d4, d5, d6]

              u44 = (1.0)*u_i[ d1, d2, d3, d4 - 2, d5, d6] +
                    (2.0)*u_i[ d1, d2, d3, d4 - 1, d5, d6] +
                    (3.0)*u_i[ d1, d2, d3, d4    , d5, d6] +
                    (4.0)*u_i[ d1, d2, d3, d4 + 1, d5, d6] +
                    (5.0)*u_i[ d1, d2, d3, d4 + 2, d5, d6]

              u55 = (1.0)*u_i[ d1, d2, d3, d4, d5 - 2, d6] +
                    (2.0)*u_i[ d1, d2, d3, d4, d5 - 1, d6] +
                    (3.0)*u_i[ d1, d2, d3, d4, d5    , d6] +
                    (4.0)*u_i[ d1, d2, d3, d4, d5 + 1, d6] +
                    (5.0)*u_i[ d1, d2, d3, d4, d5 + 2, d6]

              u66 = (1.0)*u_i[ d1, d2, d3, d4, d5, d6 - 2] +
                    (2.0)*u_i[ d1, d2, d3, d4, d5, d6 - 1] +
                    (3.0)*u_i[ d1, d2, d3, d4, d5, d6    ] +
                    (4.0)*u_i[ d1, d2, d3, d4, d5, d6 + 1] +
                    (5.0)*u_i[ d1, d2, d3, d4, d5, d6 + 2]

              u_ip1[ d1, d2, d3, d4, d5, d6] = delta_12*u11 + delta_22*u22 + delta_32*u33 + 
                                               delta_42*u44 + delta_52*u55 + delta_62*u66


            end
          end
        end
      end
    end
  end




  return nothing
end

function runtest()
  dim = 20
  arr_size = (dim^6)*8/(1024^2)
  println("array size = ", arr_size, " megabytes")

  arr1 = rand(dim, dim, dim, dim, dim, dim)
  arr2 = zeros(size(arr1))

  println("testing outer_func")
  @time outer_func(arr1, arr2)
  fill!(arr2, 0.0)
  @time outer_func(arr1, arr2)
  fill!(arr2, 0.0)

  println("testing outer_func2")
  @time outer_func2(arr1, arr2)
  fill!(arr2, 0.0)
  @time outer_func2(arr1, arr2)
  fill!(arr2, 0.0)

  dtype = typeof(arr1)
  fname = "outer_func_llvm.txt"
  f = open(fname, "w")

  println(f, "outer_func code_llvm = ")
  code_llvm(f, outer_func, (dtype, dtype))
  close(f)

  fname = "outer_func2_llvm.txt"
  f = open(fname, "w")
  println(f, "outer_func2 code_llvm = ")
  code_llvm(f, outer_func2, (dtype, dtype))
  close(f)

end

runtest()


