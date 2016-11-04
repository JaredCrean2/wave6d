function generate_setup(maxdim::Integer)

  fname = "setup2.jl"
  println("generating ", fname)
  f = open(fname, "w")

  search_func = getMpiSetup(maxdim)
  println(f, search_func)


  close(f)
  return nothing
end

function getMpiSetup(maxdim)

# comm is a MPI communicator, N is the number of dimensions


  str = "# This function performs an exhaustive search for the optimal"
  str *= "# number of processors to assign to each dimension of an N"
  str *= "# dimensional grid, where comm_size is the number of processors total"
  str *= "# Note that this procedure is too slow to use in general, but is"
  str *= "# for verification of a faster algorithm"
  str = "function getMPIMatches(comm_size::Integer, N::Integer)\n"
  str *= "# comm_size is the number of MPI ranks, N is the number of dimensions\n\n"
  indent = "  "
  str *= indent*"matches = Array(Int, 100, N)\n"
  str *= indent*"idx = 1\n"


  str *= "\n"

  for i=1:maxdim
    str *= indent*"if N == $i\n\n"
    indent  *= "  "

    search_i = getSearchN(i)
    search_i = indentString(length(indent), search_i)
    str *= search_i

    indent = indent[1:end-2]

    tmp_str = "\n"*indent*"end\n\n"
    str *= tmp_str
  end

  str *= indent*"return matches, idx - 1\n"
  str *= "end\n"

  return str
end

function getSearchN(N::Integer)

  str =""

  indent = ""
  for i=1:N
    varname = "d$i"
    str *= indent*string("for ", varname, "=1:comm_size\n")
    indent *="  "
  end

  # calculate product
  str_inner = string(indent, "val = 1")
  for i=1:N
    varname = "d$i"
    str_inner *= string("*", varname)
  end

  str *= str_inner*"\n"

  # check if value matches
  str *= "\n"*indent*"if val == comm_size\n"
  indent *= "  "

  for i=1:N
    str *= indent*"matches[idx, $i] = d$i\n"
  end
  str *= indent*"idx += 1\n"

  # check if array needs resizing
  str *= indent*"if idx > size(matches, 1)\n"
  indent *= "  "
  str *= indent*"matches = resize_arr(matches)\n"
  indent = indent[1:end-2]
  str *= indent*"end\n"


  indent = indent[1:end-2]
  str *= indent*"end\n\n"  # end if val == commsize
  indent = indent[1:end-2]

  # end all for loops
  for i=1:N
    str *= indent*"end\n"
    indent = indent[1:end-2]
  end

  return str
end

function indentString(indent::Integer, str::ASCIIString)
# indent all lines of a string a certain number of spaces
  indent_str = " "^indent
  newline_indent = "\n"*indent_str
  str = indent_str*str
  str = replace(str, "\n", newline_indent)

  return str
end


  


function getParams(maxdim::Integer)
# get the ParamType and its constructor


  str = "type ParamType{N}\n"
  indent = "  "
  
  # delta_invs, Ns, comm::MPI_Comm
  # ias, ibs
  # send, receive buffers
  # send_reqs, recv_reqs
  # xLs, xRs  -> make 2d array
  # nghost
  # coords

  # commnly used typetags
  tag_floatarr = "::Array{Float64, 1}\n"
  tag_intarr = "::Array{Int, 1}\n"

  # add fields
  str *= string(indent, "deltax_invs", tag_floatarr)
  str *= string(indent, "comm", "::MPI.Comm\n")
  str *= string(indent, "Ns", tag_intarr)
  str *= string(indent, "ias", tag_intarr)
  str *= string(indent, "ibs", tag_intarr)
  str *= string(indent, "send_reqs", "::Array{MPI.Request, 1}\n")
  str *= string(indent, "recv_reqs", "::Array{MPI.Request, 1}\n")
  str *= string(indent, "xLs", "::Array{Float64, 2}\n")
  str *= string(indent, "nghost", "::Int\n")
  str *= string(indent, "coords", "::Array{LinSpace{Float64}, 1}\n")

  str *= "end\n"

  return str
end

generate_setup(6)




