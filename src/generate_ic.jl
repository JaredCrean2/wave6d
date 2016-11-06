
function generate_ic(maxdim, prefix="")
  
  fname = prefix*"ic.jl"
  println("generating file ", fname)
  f = open(fname, "w")

  for i=1:maxdim
    str = getICFunction(i)
    println(f, str)
    println(f, "\n")
  end

  close(f)

  return nothing
end

function getICFunction(dim::Integer)

  dp1 = dim + 1
  str = ""
  str *= "function IC1{T}(params::ParamType{$dim}, u_arr::AbstractArray{T, $dp1}, t=0.0)\n"

  str *= "\n"
  indent = "  "
  # loop bounds
  for i=1:dim
    idx = dim - i + 1
    varname = string("d", i, "min")
    str *= indent*varname*" = params.ias[ $idx ]\n"

    varname = string("d", i, "max")
    str *= indent*varname*" = params.ibs[ $idx ]\n"
  end

  str *= "\n"

  # loops
  for i=1:dim
    var1 = string("d", i)
    var2 = var1*"min"
    var3 = var1*"max"
    str *= indent*"for "*var1*" = "*var2*":"*var3*"\n"
    indent *= "  "
    # get coordinate
    idx = dim - 1 + i
    str *= indent*"x$i = params.coords[ $idx ][ $var1 ]\n"
  end

  # get function
  for eq=1:2

    str_inner = ""
    # get array element
    str_inner *= getu_arr(dim, eq)
    str_inner *= " = "

    # get the terms of the function
    str_inner *= getfunc(dim, eq)

    # update string
    str *= indent*str_inner
  end


  # end statement
  for i=1:dim
    indent = indent[1:end-2]
    str *= indent*"end\n"
  end

  str *= "\n"
  str *= indent*"return nothing\n"

  indent = indent[1:end-2]

  str *= "end\n"

  return str
end

function getu_arr(dim::Integer, eq::Integer)
# get u_arr[d1, d2, ...]
# eq is the equation number

  # function evaluation
  str_inner = "u_arr[ "
  for i=1:dim
    var = string("d", i)
    str_inner *= var*", "
  end

  str_inner *= "$eq ]"

  return str_inner
end

function getfunc(dim::Integer, eq::Integer)
# ge the initial condition in dim dimensions
  str_inner = ""

  for i=1:dim
    var = "x$i"
    str_inner *= "-sin( $var + t ) + "
  end

  # remove trailing punctuation
  str_inner = str_inner[1:end-3]
  str_inner *= "\n"

  return str_inner
end
#generateIC(2)
