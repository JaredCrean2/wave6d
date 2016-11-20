type ParamType{N}
  ias::Array{Int, 1}
  ibs::Array{Int, 1}
  deltax_invs2::Array{Float64, 1}
  idxs::Array{Int, 2}
end

function runcase()

  N = 2  # dimensionality
  npoints = 64*64  # number of points
  ia = 3  # start point
  ib = npoints + ia - 1  # end point
  npoints_tot = npoints + 4  # 2 ghost points on either side

  @assert npoints == (ib - ia + 1)
  @assert npoints % 16 == 0  # npoints must be divisible by block size

  nblocks = div(npoints, 16)

  # square grid
  ias = Array(Int, N); fill!(ias, ia)
  ibs = Array(Int, N); fill!(ibs, ib)
  deltax_invs2 = Array(Float64, N); fill!(deltax_invs2, 0.5)

  # array of indices used by hilbertLoop
  idxs = Array(Int, 2, nblocks*nblocks)
  idx = 1
  for j=1:nblocks
    for i=1:nblocks
      idxs[1, idx] = i
      idxs[2, idx] = j
      idx += 1
    end
  end

  params = ParamType{N}(ias, ibs, deltax_invs2, idxs)

  u_i = rand(npoints_tot, npoints_tot, 2)
  u_ip1 = zeros(u_i)

  println("blockLoop timing:")
  @time blockLoop5_2_16(params, u_i, u_ip1)
  fill!(u_ip1, 0.0)
  @time blockLoop5_2_16(params, u_i, u_ip1)

  println("hilbertLoop timing:")
  @time hilbertLoop5_2_16(params, u_i, u_ip1)
  fill!(u_ip1, 0.0)
  @time hilbertLoop5_2_16(params, u_i, u_ip1)


  # write machine code to file
  originalstdout = STDOUT
  f = open("code_llvm_hilbert.txt", "w")
  redirect_stdout(f)

  code_llvm(hilbertLoop5_2_16, (typeof(params), typeof(u_i), typeof(u_i)))
  flush(f)
  close(f)
  
  f = open("code_native_hilbert.txt", "w")
  redirect_stdout(f)
  code_native(hilbertLoop5_2_16, (typeof(params), typeof(u_i), typeof(u_i)))

  flush(f)
  close(f)

  f = open("code_llvm_loop.txt", "w")
  redirect_stdout(f)

  code_llvm(blockLoop5_2_16, (typeof(params), typeof(u_i), typeof(u_i)))
  flush(f)
  close(f)
  
  f = open("code_native_loop.txt", "w")
  redirect_stdout(f)
  code_native(blockLoop5_2_16, (typeof(params), typeof(u_i), typeof(u_i)))

  flush(f)
  close(f)

  redirect_stdout(originalstdout)

  return nothing
end


function blockLoop5_2_16{T}(params::ParamType{2}, u_i::AbstractArray{T, 3}, u_ip1::AbstractArray{T, 3})
# decompose grid in blocks of size 16
# outer set of loops goes over blocks, inner set over point within block
# this vectorizes

  d2min = params.ias[2]
  d2max = params.ibs[2]
  blockstart2 = d2min
  nblocks2 = div( d2max - d2min+ 1, 16)
  blockend2 = blockstart2 + (16 * nblocks2) - 1
  d1min = params.ias[1]
  d1max = params.ibs[1]
  blockstart1 = d1min
  nblocks1 = div( d1max - d1min+ 1, 16)
  blockend1 = blockstart1 + (16 * nblocks1) - 1

  # only allow evenly divisible meshes
  @assert d1max == blockend1
  @assert d2max == blockend2

  @simd for d2 = blockstart2:16:blockend2
    @simd for d1 = blockstart1:16:blockend1
      @simd for d2block = 0:15
        blockidx2 = d2 + d2block
        @simd for d1block = 0:15
          blockidx1 = d1 + d1block
          idx = (blockidx1, blockidx2)
          kernel5(params, idx, u_i, u_ip1)
        end
      end
    end
  end

  return nothing
end

function hilbertLoop5_2_16{T}(params::ParamType{2}, u_i::AbstractArray{T, 3}, u_ip1::AbstractArray{T, 3})
# decompose grid into blocks of size 16
# the first loop uses params.idxs to specify the indices of the block
# the inner loops go over the points within the block

  d2offset = params.ias[2] - 1 + 1
  d1offset = params.ias[1] - 1 + 1

  idxs = params.idxs
  @simd for i=1:size(idxs, 2)
    d1 = idxs[1, i]
    d2 = idxs[2, i]

    @simd for d1block = 0:15
      d1blockidx = 16*(d1 - 1) + d1block + d1offset
      @simd for d2block = 0:15
        d2blockidx = 16*(d2 - 1) + d2block + d2offset
        idx = (d1blockidx, d2blockidx)
        kernel5(params, idx, u_i, u_ip1)
      end
    end
  end

  return nothing
end


@inline function kernel5{T}(params::ParamType{2}, idx,
                    u_i::AbstractArray{T,3}, u_ip1::AbstractArray{T,3})

  d1 = idx[1]
  d2 = idx[2]
  
  delta_12 = params.deltax_invs2[1]
  delta_22 = params.deltax_invs2[2]
  
  
  u11_1 =  (1/12)*u_i[ d1 - 2, d2, 1] +
           (-2/3)*u_i[ d1 - 1, d2, 1] +
              (0)*u_i[ d1    , d2, 1] +
            (2/3)*u_i[ d1 + 1, d2, 1] +
          (-1/12)*u_i[ d1 + 2, d2, 1]
  
  u22_1 =  (1/12)*u_i[ d1, d2 - 2, 1] +
           (-2/3)*u_i[ d1, d2 - 1, 1] +
              (0)*u_i[ d1, d2    , 1] +
            (2/3)*u_i[ d1, d2 + 1, 1] +
          (-1/12)*u_i[ d1, d2 + 2, 1]
  
  u_ip1[ d1, d2, 2] = delta_12*u11_1 + delta_22*u22_1
  
  u11_2 =  (1/12)*u_i[ d1 - 2, d2, 2] +
           (-2/3)*u_i[ d1 - 1, d2, 2] +
              (0)*u_i[ d1    , d2, 2] +
            (2/3)*u_i[ d1 + 1, d2, 2] +
          (-1/12)*u_i[ d1 + 2, d2, 2]
  
  u22_2 =  (1/12)*u_i[ d1, d2 - 2, 2] +
           (-2/3)*u_i[ d1, d2 - 1, 2] +
              (0)*u_i[ d1, d2    , 2] +
            (2/3)*u_i[ d1, d2 + 1, 2] +
          (-1/12)*u_i[ d1, d2 + 2, 2]
  
  u_ip1[ d1, d2, 1] = delta_12*u11_2 + delta_22*u22_2
  
  
  return nothing
end


runcase()
