type ParamType{N}
  delta_xs::Array{Float64, 1}
  deltax_invs2::Array{Float64, 1}
  comm::MPI.Comm
  Ns_global::Array{Int, 1}  # number of regular grid points in each direction
  Ns_local::Array{Int, 1}
  Ns_total_local::Array{Int, 1}  # number to regular grid points + ghosts in each direction
  ias::Array{Int, 1}
  ibs::Array{Int, 1}
  send_reqs::Array{MPI.Request, 1}
  recv_reqs::Array{MPI.Request, 1}
  peernums::Array{Int, 2}  # 2 x ndim array containing peer numbers
  xLs::Array{Float64, 2}  # xmin and xmax for each dimension
  nghost::Int
  coords::Array{LinSpace{Float64}, 1}
end


function ParamType(Ns_global::Array{Int, 1}, xLs::Array{Float64, 2}, nghost)
# Ns = number of grid points (not including ghosts
# xls = 2 x ndim array of xmin and xmax for each dimension

  Ns_local = copy(Ns_global)  # change when parallelizing
  Ns_total_local = Ns_local + 2*nghost

  N = length(Ns_global)
  comm = MPI.COMM_WORLD

  ias = nghost*ones(N) + 1
  ibs = zeros(ias)
  for i=1:N
    ibs[i] = Ns_total_local[i] - nghost
  end

  delta_xs = zeros(N)
  for i=1:N
    xmin = xLs[1, i]
    xmax = xLs[2, i]
    N_i = Ns[i]
    delta_xs[i] = (xmax - xmin)/(N_i - 1)
  end

  delta_xinvs2 = 1./(delta_xs.^2)

  send_reqs = Array(MPI.Request, 0)
  recv_reqs = Array(MPI.Request, 0)
  peernums = Array(Int, 0)

  # TODO update this when parallelizing
  coords = Array(LinSpace{Float64}, N)
  for i=1:N
    xmin = xls[1, i]
    xmax = xLs[2, i]
    delta_x = delta_xs[i]
    Ntot = Ns_total_local[i]

    xmin = xmin - nghost*delta_x
    xmax = xmax + nghost*delta_x

    coords[i] = linspace(xmin, xmax, Ntot)
  end

  return ParamType{N}(delta_xs, deltax_invs2, comm, Ns_global, Ns_local, Ns_total_local,
                      ias, ibs, send_reqs, recv_reqs, peernums, xLs, nghost, coords)
end







