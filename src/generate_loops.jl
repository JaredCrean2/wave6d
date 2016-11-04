# generate simple for loops for all dimensions

function generate_loops(Nmax, npts)
# Nmax is the maximum number of dimensions
# npts is the number of points in the stencil (used to figure out the name of
# the kernel to call)

  fname = string("loops_", npts, ".jl")
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

  str = ""
  str *= "function simpleLoop$npts{T}(params::ParamType{$N}, u_i::AbstractArray{T, $N}, u_ip1::AbstractArray{T, $N})\n"

  indent = "  "

  # generate loop bounds
  for i=1:N
    varname = string("d", i, "min")
    varname2 = string("d", i, "max")
    idx = N - i + 1
    str *= indent*varname*" = params.ias[$idx]\n"
    str *= indent*varname2*" = params.ibs[$idx]\n"
  end

  str *= "\n"
  # generate loops
  for i=1:N
    varname = string("d", i)
    varname_min = string("d", i, "min")
    varname_max = string("d", i, "max")

    str *= indent*"for "*varname*" = "*varname_min*":"*varname_max*"\n"
    indent *= "  "
  end

  # generate body
  str_inner = ""
  str_inner *= indent*"idx = ("
  for i=1:N
    idx = N - i + 1  # indices are reversed
    varname = "d$idx"
    str_inner *= varname*", "
  end

  # remove trailing space and comma
  str_inner = str_inner[1:end-2]
  str_inner *= ")\n"

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





generate_loops(3, 5)
