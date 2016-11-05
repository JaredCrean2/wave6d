# do communication of ghost points

function doComm{T}(params::ParamType{N}, u_arr)

  # copy data into send buffers
  for i=1:N
    for j=1:2  # upper and lower directions
      send_buf_i = params.send_bufs[j, i]

      if j == 1  # TODO: verify this
        isupper = false
      else
        isupper = true
      end

      copyToBuffer(params, u_arr, i, isupper)
   end
  end



    

