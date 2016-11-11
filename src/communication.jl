# do communication of ghost points

#=
global const TAG_SEND_UPPER = 1
global const TAG_SEND_LOWER = 2
global const TAG_RECV_LOWER = TAG_SEND_UPPER
global const TAG_RECV_UPPER = TAG_SEND_LOWER
global const SEND_TAGS = [TAG_SEND_LOWER, TAG_SEND_UPPER]
global const RECV_TAGS = [TAG_RECV_LOWER, TAG_RECV_UPPER]
=#

function startComm{N}(params::ParamType{N}, u_arr)

  # posot all receives
  for i=1:N
    for j=1:2
      peer_j = params.peernums[j, i]
      recv_buf_j = params.recv_bufs[j, i]
      recv_waited = params.recv_waited[j, i]
      recv_tag = params.recv_tags[j, i]
      recv_req = params.recv_reqs[j, i]

      if !recv_waited
        stat = MPI.Wait!(recv_req)
        if MPI.Get_error(stat) != 0 && !FORCE_SYNC
          throw(ErrorException("MPI Wait on receive $i, $j errored"))
        end
        # pedantic
        params.recv_waited[j, i] = true
      end
      recv_req = MPI.Irecv!(recv_buf_j, peer_j, recv_tag, params.comm)

      params.recv_reqs[j, i] = recv_req
      params.recv_waited[j, i] = false
    end
  end




  # copy data into send buffers and send it
  for i=1:N
#    println(params.f, "starting communication for direction ", i)
    for j=1:2  # upper and lower directions
#      println(params.f, "j = ", j)
      send_buf_j = params.send_bufs[j, i]
      peer_j = params.peernums[j, i]
      send_waited = params.send_waited[j, i]
      send_req = params.send_reqs[j, i]
      recv_req = params.recv_reqs[j, i]
      isperiodic = params.periodic_flags[j, i]
      send_tag = params.send_tags[j, i]

      if j == 1
        isupper = false
      else
        isupper = true
      end

      #TODO: optimize the case where peer_j = myrank
      copyToBuffer(params, u_arr, send_buf_j, i, isupper, isperiodic)
#=
      indices = Array(Any, N+1)
      fill!(indices, 3)
      indices[N+1] = 1
      indices[i] = 1:2

      println(params.f, "sending to process ", peer_j, ", i = ", i, ", isupper = ", isupper, ", isperiodic = ", isperiodic, ", send tag = ", send_tag)
      println(params.f, "send_buffer[1:2, 3, 3] = ", send_buf_j[indices...])
      println(params.f, "send_buffer sum = ", sum(send_buf_j))
#      println(params.f, "recv_buffer[1:2, 3, 3] = ", recv_buf_j[indices...])
      
=#
      # must wait for previous communication with peer to finish, because
      # MPI does not guarantee order of arrival
      if !send_waited
#        println(params.f, "waiting for previous send to complete")
        stat = MPI.Wait!(send_req)

        if MPI.Get_error(stat) != 0 && !FORCE_SYNC
          throw(ErrorException("MPI Wait on send $i, $j errored"))
        end

#        println(params.f, "stat error value = ", MPI.Get_error(stat))
        # pedantic
        params.send_waited[j, i] = true
      end
      send_req = MPI.Isend(send_buf_j, peer_j, send_tag, params.comm)

      params.send_reqs[j, i] = send_req
      params.send_waited[j, i] = false
    end
  end

  return nothing
end

function finishComm{N}(params::ParamType{N}, u_arr)
# wait for communication to finish and copy into u_arr

  for i=1:N
    for j=1:2
      recv_buf_j = params.recv_bufs[j, i]
      peer_j = params.peernums[j, i]
      recv_waited = params.recv_waited[j, i]
      recv_req = params.recv_reqs[j, i]

      if j == 1
        isupper = false
      else
        isupper = true
      end

      if !recv_waited
#        println(params.f, "waiting for recv ", i, ", ", j)
        stat = MPI.Wait!(recv_req)
      
        if MPI.Get_error(stat) != 0 && !FORCE_SYNC
          throw(ErrorException("MPI Wait on receive $i, $j errored"))
        end

        params.recv_waited[j, i] = true
#      else
#        println(params.f, "not waiting for recv ", i, ", ", j)
      end
#=
      indices = Array(Any, N+1)
      fill!(indices, 3)
      indices[N+1] = 1
      indices[i] = 1:2

      println(params.f, "receiving from process ", peer_j, ", recv_tag = ", params.recv_tags[j, i]) 
      println(params.f, "recv_buf[1:2, 3, 3] = ", recv_buf_j[indices...])
      println(params.f, "sum recv_buf = ", sum(recv_buf_j))
=#

      copyToMain(params, u_arr, recv_buf_j, i, isupper)
    end
  end
#=
  # for added paranoia, put bad values in ghost cells
  s1 = sub(u_arr, 43:44, :, 1:2)
  s2 = sub(u_arr, 43:44, :, 83:84)
  s3 = sub(u_arr, 43:44, 1:2, :)
  s4 = sub(u_arr, 43:44, 83:84, :)
  s5 = sub(u_arr, 1:2, :, 1:2)
  s6 = sub(u_arr, 1:2, :, 83:84)
  s7 = sub(u_arr, 1:2, 1:2, :)
  s8 = sub(u_arr, 1:2, 83:84, :)

  fill!(s1, 1e10)
  fill!(s2, 1e10)
  fill!(s3, 1e10)
  fill!(s4, 1e10)
  fill!(s5, 1e10)
  fill!(s6, 1e10)
  fill!(s7, 1e10)
  fill!(s8, 1e10)
=#
  return nothing
end


