
function generate_calcerr(maxdim, prefix="")
  
  fname = prefix*"calcerr.jl"
  println("generating file ", fname)
  f = open(fname, "w")

  for i=1:maxdim
    str = getCalcErrFunction(i)
    println(f, str)
    println(f, "\n")
  end

  close(f)

  return nothing
end

function getCalcErrFunction(dim::Integer)

  str = ""
  str *= "function calcErr1{T}(params::ParamType{$dim}, u_arr::AbstractArray{T, $dim}, t)\n"

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

  str *= indent*"max_err = 0.0\n"

  str *= "\n"

  # loops
  for i=1:dim
    var1 = string("d", i)
    var2 = var1*"min"
    var3 = var2*"max"
    str *= indent*"for "*var1*" = "*var2*":"*var3*"\n"
    indent *= "  "
    # get coordinate
    idx = dim - i + 1
    str *= indent*"x$i = params.coords[ $idx ][ $var1 ]\n"
  end

  str *= "\n"

  # function evaluation
  str_inner = "u_i = u_arr[ "
  for i=1:dim
    var = string("d", i)
    str_inner *= var*", "
  end

  # remove trailing punctuation
  str_inner = str_inner[1:end-2]
  str_inner *= "]\n"

  str_inner *= indent*"u_ex_i = "

  # the terms of the function
  for i=1:dim
    var = "x$i"
    str_inner *= "cos( $var + t ) + "
  end

  # remove trailing punctuation
  str_inner = str_inner[1:end-3]
  str_inner *= "\n"

  str_inner *= indent*"err_i = abs(u_i - u_ex_i)\n\n"
  str_inner *= indent*"if err_i > max_err\n"
  indent *= "  "
  str_inner *= indent*"max_err = err_i\n"
  indent = indent[1:end-2]
  str_inner *= indent*"end\n\n"

  str *= indent*str_inner

  # end statements
  for i=1:dim
    indent = indent[1:end-2]
    str *= indent*"end\n"
  end

  str *= "\n"
  str *= indent*"return max_err\n"

  indent = indent[1:end-2]

  str *= "end\n"

  return str
end


#generate_calcerr(2)
