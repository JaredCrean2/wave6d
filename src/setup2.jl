function getMPIMatches(comm_size::Integer, N::Integer)
# comm_size is the number of MPI ranks, N is the number of dimensions

  matches = Array(Int, 100, N)
  idx = 1

  if N == 1

    for d1=1:comm_size
      val = 1*d1
    
      if val == comm_size
        matches[idx, 1] = d1
        idx += 1
        if idx > size(matches, 1)
          matches = resize_arr(matches)
        end
      end
    
    end
    
  end

  if N == 2

    for d1=1:comm_size
      for d2=1:comm_size
        val = 1*d1*d2
    
        if val == comm_size
          matches[idx, 1] = d1
          matches[idx, 2] = d2
          idx += 1
          if idx > size(matches, 1)
            matches = resize_arr(matches)
          end
        end
    
      end
    end
    
  end

  if N == 3

    for d1=1:comm_size
      for d2=1:comm_size
        for d3=1:comm_size
          val = 1*d1*d2*d3
    
          if val == comm_size
            matches[idx, 1] = d1
            matches[idx, 2] = d2
            matches[idx, 3] = d3
            idx += 1
            if idx > size(matches, 1)
              matches = resize_arr(matches)
            end
          end
    
        end
      end
    end
    
  end

  if N == 4

    for d1=1:comm_size
      for d2=1:comm_size
        for d3=1:comm_size
          for d4=1:comm_size
            val = 1*d1*d2*d3*d4
    
            if val == comm_size
              matches[idx, 1] = d1
              matches[idx, 2] = d2
              matches[idx, 3] = d3
              matches[idx, 4] = d4
              idx += 1
              if idx > size(matches, 1)
                matches = resize_arr(matches)
              end
            end
    
          end
        end
      end
    end
    
  end

  if N == 5

    for d1=1:comm_size
      for d2=1:comm_size
        for d3=1:comm_size
          for d4=1:comm_size
            for d5=1:comm_size
              val = 1*d1*d2*d3*d4*d5
    
              if val == comm_size
                matches[idx, 1] = d1
                matches[idx, 2] = d2
                matches[idx, 3] = d3
                matches[idx, 4] = d4
                matches[idx, 5] = d5
                idx += 1
                if idx > size(matches, 1)
                  matches = resize_arr(matches)
                end
              end
    
            end
          end
        end
      end
    end
    
  end

  if N == 6

    for d1=1:comm_size
      for d2=1:comm_size
        for d3=1:comm_size
          for d4=1:comm_size
            for d5=1:comm_size
              for d6=1:comm_size
                val = 1*d1*d2*d3*d4*d5*d6
    
                if val == comm_size
                  matches[idx, 1] = d1
                  matches[idx, 2] = d2
                  matches[idx, 3] = d3
                  matches[idx, 4] = d4
                  matches[idx, 5] = d5
                  matches[idx, 6] = d6
                  idx += 1
                  if idx > size(matches, 1)
                    matches = resize_arr(matches)
                  end
                end
    
              end
            end
          end
        end
      end
    end
    
  end

  return matches, idx - 1
end

