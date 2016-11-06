# generate functions to copy data between MPI buffers (of dimensions N-1)
# and the main array (of dimension N)

function generate_buffers(nghost, Nmax, prefix="")
# nghost is the number of ghost cells on each interface
# N is the maximum number of dimenions
# prefix is prended to the file name

  fname = prefix*string("buffers_", nghost, ".jl")
  println("generating file ", fname)
  f = open(fname, "w")
  for N=1:Nmax  # loop over maximum dimensions
    str = generateFunction(nghost, N, true)
    println(f, str)
    str = generateFunction(nghost, N, false)
    println(f, str)
  end
  close(f)
  return nothing
end

function generateFunction(nghost, N, forward=true)
# generate the copy function for a particular N
# forward = true mean copy from main array to buffer
# foward = false means copy from buffer to main array

  if forward
    fname = "copyToBuffer"
  else
    fname = "copyToMain"
  end
  # write function signature
  np1 = N + 1

  if forward
    last_arg = ", isperiodic::Bool"
  else
    last_arg = ""
  end

  str = ""
  str *= "function $fname{T}(params::ParamType{$N}, u_arr::AbstractArray{T, $np1},\n           buff::AbstractArray{T, $np1}, dir::Integer, isupper::Bool $last_arg)\n"
  str *= "# copy ghost values from u_arr to buff\n"
  str *= "# dir specifies the direction to copy, (must be in the range 1:N\n"
  str *= "# isupper tells whether to copy the values for the maximum indices\n"
  str *= "# in direction dir (true), or the minimum indices\n"
  str *= "# The values are always packed in order of increasing index in \n"
  str *= "# dimension dir\n\n"

  indent = "  "
  str *= indent*"@assert params.nghost == 2\n"

  str *= "\n"


  # figure out indices of fixed indices

  # isupper case
  str *= indent*"if isupper\n"
  indent *= "  "
  ngm1 = nghost - 1
  # only the starting point differes for forward and reverse
  if forward
    str *= indent*"dfixed1 = params.ibs[dir] - $ngm1\n"
  else
    str *= indent*"dfixed1 = params.ibs[dir] + 1\n"
  end

  # rest of indices
  for i=2:nghost
    im1 = i - 1
    str *= indent*"dfixed$i = dfixed$im1 + 1\n"
  end
  indent = indent[1:end-2]

  # islower case
  str *= indent*"else\n"
  indent *= "  "

  # only the starting point differs for forward and reverse
  if forward
    str *= indent*"dfixed1 = params.ias[dir]\n"
  else
    str *= indent*"dfixed1 = 1\n"
  end
  # rest of indices
  for i=2:nghost
    im1 = i - 1
    str *= indent*"dfixed$i = dfixed$im1 + 1\n"
  end
  indent = indent[1:end-2]
  str *= indent*"end\n"


  # if this is a periodic interface, adjust indices
  # only needed for copy into buffer
  if forward 
    str *= "\n"

    # figure out the adjustment
    str *= indent*"if isperiodic\n"
    indent *= "  "
    str *= indent*"if isupper\n"
    indent *= "  "
    str *= indent*"adjust = -1\n"
    indent = indent[1:end-2]
    str *= indent*"else\n"
    indent *= "  "
    str *= indent*"adjust = 1\n"
    indent = indent[1:end-2]
    str *= indent*"end\n"
    # apply adjustment to all
    for i=1:nghost
      varname = string("dfixed", i)
      str *= indent*varname*" += adjust\n"
    end

    indent = indent[1:end-2]
    str *= indent*"end\n"
  end


  str *= "\n\n"
  for dim=1:N  # loop over excluded dimensions
    loop_indices = getindices(N, dim)

    # start if
    str *= indent*"if dir == $dim\n"
    indent *= "  "

    # figure out bounds for all loops
    # TODO: add special case for N == 1
    for j=1:(N-1)
      var_name = string("d", j, "min")
      var_name2 = string("d", j, "max")
      idx_j = loop_indices[j]
      str *= indent*var_name*" = params.ias[ $idx_j ]\n"
      str *= indent*var_name2*" = params.ibs[ $idx_j ]\n"
    end

    str *= "\n"

    # for loops over free dimensions
    for j=1:(N-1)
      str *= indent*"for d$j = "*string("d", j, "min:", "d", j, "max")*"\n"
      indent *= "  "
    end

    for j=1:nghost
      if forward
        str *= indent*getAssignment(N, dim, "dfixed$j", j, 1)

        str *= indent*getAssignment(N, dim, "dfixed$j", j, 2)
      else
        str *= indent*getAssignment2(N, dim, "dfixed$j", j, 1)
        str *= indent*getAssignment2(N, dim, "dfixed$j", j, 2)
      end
    end

    # end statements
    for j=1:(N-1)
      indent = indent[1:end-2]
      str *= indent*"end\n"
    end

    str *= indent*"\n"
    indent = indent[1:end-2]
    str *= indent*"end  #end if \n\n"


  end

  str *= indent*"return nothing\n"
  indent = indent[1:end-2]
  str *= "end\n\n"

  return str

end

function getindices(N, dir)
# get the dimensions of the dimensions that are being looped over
# N is the number of dimensions, dir is the one to exclude

  vals = collect(1:N)
  vals[dir] = 0
  sort!(vals, rev=true)  # loop over last index first
  vals_extract = vals[1:end-1]

  return vals_extract
end

function getAssignment(N, dir::Integer, dir_var::ASCIIString, ind::Integer, eq::Integer)
# assign an element of u_arr to buff
# dir_var is the name of the fixed variable in u_arr
# ind is the fixed index in buff
# eq is the equation number

  str = ""

  str *= "buff[ "
  curr_idx = 1
  for i=1:N
    if i == dir
      str *= string(ind)*", "
    else
      tmp = N - curr_idx
      str *= "d$tmp"*", "
      curr_idx += 1
    end
  end

  str *= "$eq ] = "

  str *= "u_arr[ "

  curr_idx = 1
  for i=1:N
    if i == dir
      str *= dir_var*", "
    else
      tmp = N - curr_idx
      var_name = "d$tmp"
      curr_idx += 1
      str *= var_name*", "
    end
  end

  str *= "$eq ]\n"

  return str
end

function getAssignment2(N, dir::Integer, dir_var::ASCIIString, ind::Integer, eq::Integer)
# assign an element of u_arr to buff
# dir_var is the name of the fixed variable in u_arr
# ind is the fixed index in buff
# eq is the equation number

  str = ""

  str *= "u_arr[ "

  curr_idx = 1
  for i=1:N
    if i == dir
      str *= dir_var*", "
    else
      tmp = N - curr_idx
      var_name = "d$tmp"
      curr_idx += 1
      str *= var_name*", "
    end
  end

  str *= "$eq ] = "


  str *= "buff[ "
  curr_idx = 1
  for i=1:N
    if i == dir
      str *= string(ind)*", "
    else
      tmp = N - curr_idx
      str *= "d$tmp"*", "
      curr_idx += 1
    end
  end


  str *= "$eq ]\n"

  return str
end

#generate_buffers(2, 6)
