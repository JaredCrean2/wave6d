# generate functions to copy data between MPI buffers (of dimensions N-1)
# and the main array (of dimension N)

function generateBuffers(nghost, Nmax)
# nghost is the number of ghost cells on each interface
# N is the maximum number of dimenions

  for N=1:Nmax  # loop over maximum dimensions
    for dim=1:N # loop over interfaces on this dimension


  return nothing
end

function generateFunction(nghost, N)
# generate the copy function for a particular N

  #TODO: add an if statememnt to control copy direction
  #      ie. buffer to main or main to buffer
  # write function signature
  str = ""
  str *= "function copyToBuffer{N}(params::ParamType{T,N}, u_arr::AbstractArray{T, N, N2}, buff::AbstractArray{T, N2}, dir::Integer, isupper::Bool)\n"
  str *= "# copy ghost values from u_arr to buff\n"
  str *= "# dir specifies the direction to copy, (must be in the range 1:N\n"
  str *= "# isupper tells whether to copy the values for the maximum indices\n"
  str *= "# in direction dir (true), or the minimum indices\n"
  str *= "# The values are always packed in order of increasing index in \n"
  str *= "# dimensions dir\n\n"

  for dim=1:N  # loop over interface dimensions
    # start if

    # write loops

    # end if
  end

  return str

end
