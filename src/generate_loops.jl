# generate simple for loops for all dimensions

function generate_loops(Nmax, npts, prefix="")
# Nmax is the maximum number of dimensions
# npts is the number of points in the stencil (used to figure out the name of
# the kernel to call)
# prefix is prepended to the file name

  fname = prefix*string("loops_", npts, ".jl")
  println("generating file ", fname)
  f = open(fname, "w")
  for i=1:Nmax
    str = generateLoops(i, npts)
    println(f, str)
  end

  close(f)

  return nothing
end


function generateLoops(N, npts)
# generate a function that does N nested for loops

  np1 = N + 1
  str = ""
  str *= "function simpleLoop$npts{T}(params::ParamType{$N}, u_i::AbstractArray{T, $np1}, u_ip1::AbstractArray{T, $np1})\n"

  str *= "\n"
  indent = "  "

  # generate loop bounds
  for i=N:-1:1
    varname = string("d", i, "min")
    varname2 = string("d", i, "max")
    str *= indent*varname*" = params.ias[$i]\n"
    str *= indent*varname2*" = params.ibs[$i]\n"
  end

  str *= "\n"
  # generate loops
  for i=N:-1:1
    varname = string("d", i)
    varname_min = string("d", i, "min")
    varname_max = string("d", i, "max")

    str *= indent*"@simd for "*varname*" = "*varname_min*":"*varname_max*"\n"
    indent *= "  "
  end

  str_inner, indent = get_tuple(indent, N)
  str *= str_inner

  # call kernel
  str *= indent*"kernel$npts(params, idx, u_i, u_ip1)\n"

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

function get_tuple(indent, N)
  # generate body
  str_inner = ""
  str_inner *= indent*"idx = ("
  for i=1:N
    varname = "d$i"
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
