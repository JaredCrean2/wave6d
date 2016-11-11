
# LSERK coefficients
global const a_coeffs = [0; -567301805773.0/1357537059087.0; -2404267990393.0/2016746695238.0; -3550918686646.0/2091501179385.0; -1275806237668.0/842570457699.0]

global const b_coeffs = [1432997174477.0/9575080441755.0; 5161836677717.0/13612068292357.0; 1720146321549.0/2090206949498.0; 3134564353537.0/4481467310338.0; 2277821191437.0/14882151754819.0]

function lserk(f::Function, tmax, u::AbstractArray, params)
# 4th order classical Runga-Kutta time stepping
# u is ncomp x params.Ntot matrix, where ncomp is 
# the number of components of the system of equations

# the initial conditions for the system is stored in u, and the final 
# solution is returned in u

  delta_t = params.delta_t
  tfac = delta_t/2
  tfac6 = delta_t/6
  t = 0.0

  # allocate storage
  # this is actually a 3N scheme because the function f cannot overwrite its
  # input, it needs a separate vector for output
  dU = zeros(u)
  F_vals = zeros(u)
#  u_tmp = zeros(u)
#  k2 = zeros(u)
#  k3 = zeros(u)
#  k4 = zeros(u)
#  u_tmp = zeros(u)

  tstep = 1
  while (t < tmax)

#    if tstep % 100 == 0
     if params.comm_rank == 0
       println("t = ", t)
     end
#   end
#    err = calcError(params, u, t)
#    println("  err = ", err)


    # stage 1
    # TODO: get c_i coefficients to update t
    t_i = t
    f(params, u, F_vals, t_i)

    fac = b_coeffs[1]
    for j=1:length(u)
      dU[j] = delta_t*F_vals[j]
      u[j] += fac*dU[j]
    end

    for i=2:5
      f(params, u, F_vals, t_i)
      # update
      fac = a_coeffs[i]
      fac2 = b_coeffs[i]
      for j=1:length(u)
        dU[j] = fac*dU[j] + delta_t*F_vals[j]
        u[j] += fac2*dU[j]
      end
    end

    # update
    t += delta_t
    tstep += 1

  end  # end while loop

  return t

end  # end function

