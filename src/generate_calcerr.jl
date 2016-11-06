
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

  dp1 = dim + 1
  str = ""
  str *= "function calcErr1{T}(params::ParamType{$dim}, u_arr::AbstractArray{T, $dp1}, t)\n"

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
    var3 = var1*"max"
    str *= indent*"for "*var1*" = "*var2*":"*var3*"\n"
    indent *= "  "
    # get coordinate
    idx = dim - i + 1
    str *= indent*"x$i = params.coords[ $idx ][ $var1 ]\n"
  end

  str *= "\n"

  # function evaluation
  for eq = 1:2
    str_inner = "u_i_$eq = "
    str_inner *= getu_arr(dim, eq)*"\n"

    str_inner *= indent*"u_ex_i_$eq = "
    str_inner *= getfunc(dim, eq)


    str_inner *= indent*"err_i_$eq = abs(u_i_$eq - u_ex_i_$eq)\n\n"
    str_inner *= indent*"if err_i_$eq > max_err\n"
    indent *= "  "
    str_inner *= indent*"max_err = err_i_$eq\n"
    indent = indent[1:end-2]
    str_inner *= indent*"end\n\n"

    str *= indent*str_inner
  end

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
