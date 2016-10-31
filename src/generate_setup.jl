function generate_setup(maxdim::Integer)

  fname = "setup2.jl"
  println("generating ", fname)
  f = open(fname, "w")

  params_def = getParams(maxdim)
  println(f, params_def)


  close(f)
  return nothing
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




