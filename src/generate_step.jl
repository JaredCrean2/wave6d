function generate_step(Nmax, npts, blocksizes::AbstractArray, prefix="")
# Nmax is the maximum number of dimensions
# npts is the number of points in the stencil (used to figure out the name of
# the kernel to call)
# prefix is prepended to the file name

  fname = prefix*string("step_", npts, ".jl")
  println("generating file ", fname)
  f = open(fname, "w")

  println(f, "# naming convention: number of points, number of blocked loops, blocksize")

  # do the regular blocked loops
  fname_outer = "step"
  fname_inner = "blockLoop"
#  for i=1:Nmax  # dimensions
    for j=1:Nmax  # blocked loops
      for blocksize in blocksizes
         str = generateStep(j, blocksize, npts, fname_outer, fname_inner)
         println(f, str)
       end
     end
#  end

  # now do the hilbert blocked loops
  fname_outer = "stepH"
  fname_inner = "hilbertLoop"
  for i=1:Nmax
    for blocksize in blocksizes
      str = generateStep(i, blocksize, npts, fname_outer, fname_inner)
      println(f, str)
    end
  end

   str = generate_accessor(Nmax, blocksizes, npts)
   println(f, str)

  close(f)

  return nothing
end

function generateStep(Nblock, blocksize, npts, fname_outer, fname_inner)


  str = ""
  indent = ""

  # naming convention: number of points, number of blocked loops, blocksize
  suffix = string(npts, "_", Nblock, "_", blocksize)
  str *= "function $fname_outer$suffix{N}(params::ParamType{N}, u_i, u_ip1, t)\n"
  str *= "# single timestep\n"
  indent *= "  "
  str *= indent*"params.t = t\n"

  str *= indent*"params.time.t_comm += @elapsed begin\n"
  indent *= "  "

  str *= indent*"startComm(params, u_i)\n"

  str *= indent*"finishComm(params, u_i)\n"

  str *= indent*"if FORCE_SYNC\n"
  indent *= "  "
  str *= indent*"MPI.Barrier(params.comm)\n"
  indent = indent[1:end-2]
  str *= indent*"end\n"
  indent = indent[1:end-2]

  str *= indent*"end\n"


  str *= indent*"params.time.t_compute += @elapsed $fname_inner$suffix(params, u_i, u_ip1)\n"

  str *= indent*"params.itr += 1\n"

  str *= indent*"return nothing\n"
  indent = indent[1:end-2]
  str *= indent*"end\n"

  return str
end



function generate_accessor(Nblock, blocksizes::AbstractArray, npts)

  # N doesn't matter, it is handled by dispatch
  str = ""
  indent = ""

  str *= "function getStepFunc(N::Integer, Nblock::Integer, blocksize::Integer, npts::Integer)\n"

  indent *= "  "

  # special case: no blocking
  str *= indent*"if Nblock == 0\n"
  indent *= "  "
  str *= indent*"return step, simpleLoop5\n"  # return the simple loop version
  indent = indent[1:end-2]

  # special case: no block hilbert
  str *= indent*"elseif Nblock == -1 && blocksize == 0\n"
  indent *= "  "
  str *= indent*"return hilbertStep, hilbertLoop5\n"
  indent = indent[1:end-2]

  # regular block loops
  fname_outer = "step"
  fname_inner = "blockLoop"
  for Nblock_i = 1:Nblock  # loop over number of blocked loops
    str *= indent*"elseif Nblock == $Nblock_i \n\n"
    indent *= "  "
    str *= get_blocksize_loop(indent, Nblock_i, blocksizes, npts, fname_outer, fname_inner)
    indent = indent[1:end-2]
  end

  # hilbert block loops
  fname_outer = "stepH"
  fname_inner = "hilbertLoop"
  for Nblock_i = 1:Nblock  # should be N really, because we don't support
                           # partial loop blocking for Hilbert
    str *= indent*"elseif Nblock == -1 && N == $Nblock_i\n"
    indent *= "  "
    str *= get_blocksize_loop(indent, Nblock_i, blocksizes, npts, fname_outer, fname_inner)
    indent = indent[1:end-2]
  end

  str *= indent*"else\n"
  indent *= "  "
  str *= indent*"throw(ErrorException(\"unsupported mode: Nblock = \$Nblock, blocksize = \$blocksize\"))\n"

  indent = indent[1:end-2]
  str *= indent*"end\n"

  str *= indent*"return nothing\n"
  indent = indent[1:end-2]
  str *= indent*"end\n"

  return str
end

function get_blocksize_loop(indent, Nblock, blocksizes::AbstractArray, npts, fname_outer, fname_inner)

  str = ""
  for i=1:length(blocksizes)
    blocksize_i = blocksizes[i]
    if i == 1
      if_type = "if "
    else
      if_type = "elseif "
    end

    str *= indent*if_type*"blocksize == $blocksize_i\n"
    indent *= "  "
    fname = string(fname_outer, npts, "_", Nblock, "_", blocksize_i)
    fname2 = string(fname_inner, npts, "_", Nblock, "_", blocksize_i)
    str *= indent*"return "*fname*", "*fname2*"\n"
    indent = indent[1:end-2]
  end

  str *= indent*"else\n"
  indent *= "  "
  str *= indent*"throw(ErrorException(\"unsupported mode: Nblock = \$Nblock, blocksize = \$blocksize\"))\n"

  indent = indent[1:end-2]
  str *= indent*"end\n\n"

  return str
end





