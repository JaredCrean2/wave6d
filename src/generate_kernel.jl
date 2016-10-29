# generate the spatial kernel for an n-dimensional wave equation solver

function generate_kernel(maxdim, stencil)
# writes to file kernel_dim_order

  npts = length(stencil)
  order = npts - 1
  fname = string("kernel_", maxdim, "_", order, ".jl")
  println("creating file ", fname)
  f = open(fname, "w")

  second_line_indent = " "^6
  for dim=1:maxdim
    var_name = string("u", dim, dim)
    stencil_dim = getStencil(maxdim, stencil, dim, second_line_indent)
    println(f, var_name, " = ", stencil_dim)
  end

  kernel_assembly = getKernelAssembly(maxdim)
  println(f, kernel_assembly)

  close(f)

end


function getStencil(maxdim::Integer, stencil, dim::Integer, 
                    second_line_indent::ASCIIString="")

# Generates a string that multiplies the given stencil coefficients by
#   the corresponding values from an array called u_i.
# The stencil coefficients are list from lowest index to highest index
#   maxdim is the number of indices the array u_i has
# The indices are named d1 through dn.
# dim is the dimension of u_i upon which to apply the array.
# Each term of the stencil is put on its own line.  second_line_indent
# is the string prepended to all lines except the first to indent them

# the stencil array is expected to be an array of strings
# the strings are enclosed in parenthesis before multiplication
# great pains are taken to ensure alignment of all terms of the stencil
  
  npts = length(stencil)
  stencil_startoffset = -div(npts - 1, 2)
  stencil_endoffset = div(npts - 1, 2)
  stencil_offsets = stencil_startoffset:stencil_endoffset

  # figure out maximum size of a stencil coefficient
  max_coeff_len = 0
  for i=1:npts
    len_i = length(stencil[i])
    if len_i > max_coeff_len
      max_coeff_len = len_i
    end
  end


  str = ""

  # index names are d1 through dmaxdim
  for pt=1:npts

    if pt == 1
      indentlevel = ""
    else
      indentlevel = second_line_indent
    end

    stencil_pt = stencil[pt]
    stencil_pt_len = length(stencil_pt)
    stencil_offset = stencil_offsets[pt]

    padding = max_coeff_len - stencil_pt_len
    coeff_str = string( " "^padding, "(", stencil_pt, ")")

    str_pt = indentlevel*coeff_str*"*u_i[ "
    
    # insert all indices before dim
    for i=1:(dim-1)
      idx_i = "d$i"
      str_pt = str_pt*idx_i*", "
    end

    # insert index dim and modifier
    
    # figure out if this should be a plus or minus
    stencil_offset_inner = abs(stencil_offset)
    if stencil_offset < 0
      mod_str = "d$dim - $stencil_offset_inner, "
    elseif stencil_offset > 0
      mod_str = "d$dim + $stencil_offset_inner, "
    else  # stencil_offset == 0
      mod_str = "d$dim    , "
    end

    # apply the index dim and modifier
    str_pt *= mod_str


    # insert all indices after dim
    for i=(dim+1):maxdim
      idx_i = "d$i"
      str_pt *= idx_i*", "
    end

    # remove last comma and space
    str_pt = str_pt[1:end-2]

    # close bracket
    str_pt = str_pt*"]"


    # figure out whether or not to add a plus sign for another term
    if pt != npts
      mod_str = " +\n"
    else
      mod_str = "\n"
    end

    str *= str_pt*mod_str

  end  # end loop over npts

  return str

end

function getKernelAssembly(maxdim::Integer)

  second_line_indent = " "^(9 + 4*maxdim)
  # write receiving location
  str = "u_ip1[ "
  for i=1:maxdim
    str *= "d$i, "
  end

  # remove final comma and space, add close brack
  str = str[1:end-2]
  str *= "] = "

  # multiply stencil in each direction by grid spacing square
  for i=1:maxdim
    delta_name = string("delta_", i, "2")
    stencil_fac_name = string("u", i, i)
    str *= delta_name*"*"*stencil_fac_name*" + "

    if ((i % 3) == 0) && i != maxdim
      str *= "\n"*second_line_indent
    end

  end

  # remove last " + "
  str = str[1:end-3]


  return str
end

function indentString(indent::Integer, str::ASCIIString)

  indent_str = " "^indent
  newline_indent = "\n"*indent_str
  str = indent_str*str
  str = replace(str, "\n", newline_indent)

  return str
end

stencil = ["1.0", "2", "3", "4", "5"]
maxdim = 6

generate_kernel(maxdim, stencil)
#str = getStencil(maxdim, stencil, 2)
#println("str = \n", str)
