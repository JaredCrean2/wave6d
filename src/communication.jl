# do communication of ghost points

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
    for j=1:2  # upper and lower directions
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
      
      # must wait for previous communication with peer to finish, because
      # MPI does not guarantee order of arrival
      if !send_waited
        stat = MPI.Wait!(send_req)

        if MPI.Get_error(stat) != 0 && !FORCE_SYNC
          throw(ErrorException("MPI Wait on send $i, $j errored"))
        end

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
        stat = MPI.Wait!(recv_req)
      
        if MPI.Get_error(stat) != 0 && !FORCE_SYNC
          throw(ErrorException("MPI Wait on receive $i, $j errored"))
        end

        params.recv_waited[j, i] = true
      end
      
      copyToMain(params, u_arr, recv_buf_j, i, isupper)
    end
  end

  return nothing
end


