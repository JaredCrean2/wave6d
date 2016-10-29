# test whether it is necessary to manually inline the kernel into to function
# containing the loops
# this is a testcase for the indexing slowdown in julia 0.5/0.6

function outer_func2(u_i, u_ip1)
  s1 = size(u_i, 1)
  nghost = 2
  ia = nghost + 1
  ib = s1 - nghost

  for d6=ia:ib
    for d5=ia:ib
      for d4=ia:ib
        for d3 = ia:ib
          for d2= ia:ib
            @simd for d1 = ia:ib


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
  # dim = 20 is used for benchmarking
  dim = 15
  arr_size = (dim^6)*8/(1024^2)
  println("array size = ", arr_size, " megabytes")

  arr1 = rand(dim, dim, dim, dim, dim, dim)
  arr2 = zeros(size(arr1))

  println("testing outer_func2")
  f= open("profileout.dat", "w")
  @time outer_func2(arr1, arr2)
  fill!(arr2, 0.0)
  @profile @time outer_func2(arr1, arr2)
  Profile.print(f, format=:flat)
  fill!(arr2, 0.0)

#=
  dtype = typeof(arr1)

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


