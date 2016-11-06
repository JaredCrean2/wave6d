# do communication of ghost points

global const TAG_SEND_UPPER = 1
global const TAG_SEND_LOWER = 2
global const TAG_RECV_LOWER = TAG_SEND_UPPER
global const TAG_RECV_UPPER = TAG_SEND_LOWER
global const SEND_TAGS = [TAG_SEND_LOWER, TAG_SEND_UPPER]
global const RECV_TAGS = [TAG_RECV_LOWER, TAG_RECV_UPPER]

function startComm{N}(params::ParamType{N}, u_arr)

  # copy data into send buffers
  for i=1:N
    for j=1:2  # upper and lower directions
      send_buf_j = params.send_bufs[j, i]
      recv_buf_j = params.recv_bufs[j, i]
      peer_j = params.peernums[j, i]
      send_waited = params.send_waited[j, i]
      recv_waited = params.recv_waited[j, i]
      send_req = params.send_reqs[j, i]
      recv_req = params.recv_reqs[j, i]

      if j == 1
        isupper = false
      else
        isupper = true
      end

      #TODO: optimize the case where peer_j = myrank
      copyToBuffer(params, u_arr, send_buf_j, i, isupper)

      # must wait for previous communication with peer to finish, because
      # MPI does not guarantee order of arrival
      if !send_waited
        MPI.Wait!(send_req)
        # pedantic
        params.send_waited[j, i] = true
      end
      send_req = MPI.Isend(send_buf_j, peer_j, SEND_TAGS[j], params.comm)

      if !recv_waited
        MPI.Wait!(recv_req)
        # pedantic
        params.recv_waited[j, i] = false
      end
      recv_req = MPI.Irecv!(recv_buf_j, peer_j, RECV_TAGS[j], params.comm)

      params.send_reqs[j, i] = send_req
      params.send_waited[j, i] = false
      params.recv_reqs[j, i] = recv_req
      params.recv_waited[j, i] = false

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
        MPI.Wait!(recv_req)
        params.recv_waited[j, i] = true
      end

      copyToMain(params, u_arr, recv_buf_j, i, isupper)
    end
  end

  return nothing
end









    

