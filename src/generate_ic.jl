
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

  str = ""
  str *= "function IC1{T}(params::ParamType{$dim}, u_arr::AbstractArray{T, $dim})\n"

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

  # function evaluation
  str_inner = "u_arr[ "
  for i=1:dim
    var = string("d", i)
    str_inner *= var*", "
  end

  # remove trailing punctuation
  str_inner = str_inner[1:end-2]
  str_inner *= "]"

  str_inner *= " = "

  # the terms of the function
  for i=1:dim
    var = "x$i"
    str_inner *= "cos( $var ) + "
  end

  # remove trailing punctuation
  str_inner = str_inner[1:end-3]
  str_inner *= "\n"

  str *= indent*str_inner

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



#generateIC(2)
