# generate simple for loops for all dimensions

function generate_loops_hilbert(Nmax, npts, prefix="")
# Nmax is the maximum number of dimensions
# npts is the number of points in the stencil (used to figure out the name of
# the kernel to call)
# prefix is prepended to the file name

  fname = prefix*string("hilbertloops_", npts, ".jl")
  println("generating file ", fname)
  f = open(fname, "w")
  for i=1:Nmax
    str = generateLoopsHilbert(i, npts)
    println(f, str)
  end

  close(f)

  return nothing
end


function generateLoopsHilbert(N, npts)
# generate a function that does N nested for loops

  np1 = N + 1
  str = ""
  str *= "function hilbertLoop$npts{T}(params::ParamType{$N}, u_i::AbstractArray{T, $np1}, u_ip1::AbstractArray{T, $np1})\n"

  str *= "\n"
  indent = "  "

  # generate loop bounds
  for i=N:-1:1
    varname = string("d", i, "offset")
#    varname2 = string("d", i, "max")
    str *= indent*varname*" = params.ias[$i] - 1\n"
#    str *= indent*varname2*" = params.ibs[$i]\n"
  end

  str *= "\n"

  str *= indent*"idxs = params.idxs\n"

  str *= indent*"$macro_name for i=1:size(idxs, 2)\n"
  indent *= "  "

  str_inner, indent = get_tuple_hilbert(indent, N)
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

function get_tuple_hilbert(indent, N)
  # generate body
  str_inner = ""
  str_inner *= indent*"idx = ("
  for i=1:N
    offset_name = string("d", i, "offset")
    varname = "idxs[ $i, i] + $offset_name"
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
