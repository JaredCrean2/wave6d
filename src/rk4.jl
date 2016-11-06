
function rk4(f::Function, tmax, u::AbstractArray, params)
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

  k1 = zeros(u)
  k2 = zeros(u)
  k3 = zeros(u)
  k4 = zeros(u)
  u_tmp = zeros(u)

  tstep = 1
  while (t < tmax)

#    if tstep % 100 == 0
     if params.comm_rank == 0
       println("t = ", t)
     end
#   end
#    err = calcError(params, u, t)
#    println("  err = ", err)

    # copy u into u_tmp
    for i=1:length(u)
      u_tmp[i] = u[i]
    end

    # stage 1
    t_i = t
    f(params, u_tmp, k1, t_i)  # put result in t

    # stage 2
    t_i = t + tfac
    # calculate state for stage 2
    for i = 1:length(u)
      u_tmp[i] = u[i] + tfac*k1[i]
    end
    f(params, u_tmp, k2, t_i)

    # stage 3
    t_i = t + tfac
    for i=1:length(u)
      u_tmp[i] = u[i] + tfac*k2[i]
    end
    f(params, u_tmp, k3, t_i)

    # stage 4
    t_i = t + delta_t
    for i=1:length(u)
      u_tmp[i] = u[i] + delta_t*k3[i]
    end
    f(params, u_tmp, k4, t_i)

    # update
#    println("k1 = \n", k1)
#    println("k2 = \n", k2)
#    println("k3 = \n", k3)
#    println("k4 = \n", k4)
    t += delta_t
    for i=1:length(u)
      u[i] = u[i] + tfac6*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i])
    end

    tstep += 1

  end  # end while loop

  return t

end  # end function


type ParamsRK4
  delta_t::Float64
end

function testrk4()

  params = ParamsRK4(0.1)
  u = zeros(1, 1)
  u[1] = u_exact(0.0, 1.0)

  tmax = 1.0
  tfinal = rk4(testfunc, tmax, u, params)

  println("u = ", u[1])
  @assert abs(u[1] - u_exact(tfinal)) < 1e-12
end

function testfunc(params::ParamsRK4, u_i, t, u_ip1)
  u_ip1[1] = inner_func(t)
end

# derivataive of function
function inner_func(t)
  t_pert = t + 1e-20im
  return imag(u_exact(t_pert))/1e-20
end

# function itself
function u_exact(t)
  return t + 1
end

#testrk4()
