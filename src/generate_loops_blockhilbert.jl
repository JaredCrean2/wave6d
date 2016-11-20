# generate simple for loops for all dimensions

function generate_loops_blockhilbert(Nmax, npts, blocksize, prefix="")
# Nmax is the maximum number of dimensions
# npts is the number of points in the stencil (used to figure out the name of
# the kernel to call)
# prefix is prepended to the file name

  fname = prefix*string("hilbertblockloops_", npts, "_", blocksize, ".jl")
  println("generating file ", fname)
  f = open(fname, "w")
  for i=1:Nmax
    str = generateLoopsBlockHilbert(i, npts, blocksize)
    println(f, str)
  end

  close(f)

  return nothing
end


function generateLoopsBlockHilbert(N, npts, blocksize)
# generate a function uses the hilbert indices as the indices of blocks of 
# points of size blocksize

  np1 = N + 1
  str = ""
  # for consistency, include the number of blocked loops in the function name
  fname = string("hilbertLoop", "_", npts, "_", N, "_", blocksize)
  str *= "function $fname{T}(params::ParamType{$N}, u_i::AbstractArray{T, $np1}, u_ip1::AbstractArray{T, $np1})\n"

  str *= "\n"
  indent = "  "

  # generate loop bounds
  for i=N:-1:1
    varname = string("d", i, "offset")
#    varname2 = string("d", i, "max")
    str *= indent*varname*" = params.ias[$i] - 1 + 1\n"
#    str *= indent*varname2*" = params.ibs[$i]\n"
  end

  str *= "\n"

  str *= indent*"idxs = params.idxs\n"

  str *= indent*"@simd for i=1:size(idxs, 2)\n"
  indent *= "  "

  # unpack the indices of the bottom left corner of the current block
  for i=1:N
    varname = string("d", i)
    str *= indent*varname*" = idxs[$i, i]\n"
  end

  str *= "\n"

  # add loops over the block and calculate the true index
  for i=1:N
    outer_idxname = string("d", i)
    idxname = string("d", i, "block")  # block loop index
    rng = string("0:", blocksize - 1)
    varname = string("d", i, "blockidx")  # true index
    offset_name = string("d", i, "offset")
    str *= indent*"for $idxname = $rng\n"
    indent *= "  "
    str *= indent*varname*" = $blocksize*($outer_idxname - 1) + $idxname + $offset_name\n"
  end


  str_inner, indent = get_tuple_blockhilbert(indent, N)
  str *= str_inner

  # call kernel
  str *= indent*"kernel$npts(params, idx, u_i, u_ip1)\n"

  indent = indent[1:end-2]
  str *= indent*"end\n"

  str *= "\n"
  str *= indent*"return nothing\n"
  indent = indent[1:end-2]
  str *= "end\n"

  return str
end

function get_tuple_blockhilbert(indent, N)
  # generate body
  str_inner = ""
  str_inner *= indent*"idx = ("
  for i=1:N
    varname = string("d", i, "blockidx")
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
