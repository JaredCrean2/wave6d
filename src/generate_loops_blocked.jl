# generate blocked for loops for all dimensions

function generate_loops_blocked(Nmax, npts, blocksize, prefix="")
# Nmax is the maximum number of dimensions
# npts is the number of points in the stencil (used to figure out the name of
# the kernel to call)
# prefix is prepended to the file name

  fname = prefix*string("blockloops_", npts, "_", blocksize, ".jl")
  println("generating file ", fname)
  f = open(fname, "w")
  for i=1:Nmax
    for j=1:i  # how many loops to block
      str = generateLoops_blocked(i, npts, j, blocksize)
      println(f, str)
    end
  end

  close(f)

  return nothing
end


function generateLoops_blocked(N, npts, Nblock, blocksize)
# generate a function that does N nested for loops
# block loops 1:Nblock
# blocksize is the blocksize

  @assert Nblock <= N

  np1 = N + 1
  str = ""

  # naming convention: number of points, nuber of blocked loops, blocksize
  fname = string("blockLoop", npts, "_", Nblock, "_", blocksize)
  str *= "function $fname{T}(params::ParamType{$N}, u_i::AbstractArray{T, $np1}, u_ip1::AbstractArray{T, $np1})\n"

  str *= "\n"
  indent = "  "

  # generate loop bounds
  for i=N:-1:1
    varname = string("d", i, "min")
    varname2 = string("d", i, "max")
    str *= indent*varname*" = params.ias[$i]\n"
    str *= indent*varname2*" = params.ibs[$i]\n"

    if i <= Nblock  # do blocking
      block_start = string("blockstart", i)
      nblocks = string("nblocks", i)
      block_end = string("blockend", i)

      block_start_calc = " = "*varname*"\n"
      nblocks_calc = string(" = div( ", varname2, " - ", varname, "+ 1, ", blocksize, ")\n")
      block_end_calc = " = "*block_start*" + ("*string(blocksize)*" * "*nblocks*") - 1\n"

      str *= indent*block_start*block_start_calc
      str *= indent*nblocks*nblocks_calc
      str *= indent*block_end*block_end_calc
    end
  end

  str *= "\n"

  str *= indent*"# only allow evenly divisible meshes\n"
  for i=1:Nblock
    varname = string("d", i, "max")
    varname2 = string("blockend", i)
    str *= indent*"@assert "*varname*" == "*varname2*"\n"
  end

  str *= "\n"
  # generate loops
  for i=N:-1:1
    varname = string("d", i)
    if i <= Nblock
      varname_min = string("blockstart", i)
      varname_max = string("blockend", i)
      loop_range = string(varname_min, ":", blocksize, ":", varname_max)
    else
      varname_min = string("d", i, "min")
      varname_max = string("d", i, "max")
      loop_range = string(varname_min, ":", varname_max)
    end

    str *= indent*"@simd for "*varname*" = "*loop_range*"\n"
    indent *= "  "
  end


  # add block loops
  for i=Nblock:-1:1
    varname = string("d", i, "block")
    varname_outer = string("d", i)
    blockidx = string("blockidx", i)
    loop_range = string(0, ":", blocksize - 1)

    str *= indent*"@simd for "*varname*" = "*loop_range*"\n"
    indent *= "  "
    str *= indent*blockidx*" = "*varname_outer*" + "*varname*"\n"
  end

  str_inner, indent = get_tuple(indent, N, Nblock)
  str *= str_inner

  # call kernel
  str *= indent*"kernel$npts(params, idx, u_i, u_ip1)\n"

  # end block loops
  for i=Nblock:-1:1
    indent = indent[1:end-2]
    str *= indent*"end\n"
  end


  # end loops
  for i=1:N
    indent = indent[1:end-2]
    str *= indent*"end\n"
  end

  str *= "\n"
  str *= indent*"return nothing\n"
  indent = indent[1:end-2]
  str *= "end\n"

  return str
end

function get_tuple(indent, N, Nblock)
  # generate body
  str_inner = ""
  str_inner *= indent*"idx = ("
  for i=1:N
    if i <= Nblock
      varname = "blockidx$i"
    else
      varname = "d$i"
    end
    str_inner *= varname*", "
  end

  # remove trailing space and comma
  if N > 1
    str_inner = str_inner[1:end-2]
  end
  str_inner *= ")\n"

  return str_inner, indent
end




#generate_loops(6, 5)
