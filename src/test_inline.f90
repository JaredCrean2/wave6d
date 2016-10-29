PROGRAM test_inline
#ifdef __INTEL_COMPILER
    USE IFPORT
#endif
IMPLICIT NONE

  integer i;
  REAL*8 id;
  REAL*8 d2;
  integer, PARAMETER :: sz = 20;
  REAL*8 tstart, tfinish
  REAL*8 arr1(sz, sz, sz, sz, sz, sz);
  REAL*8 arr2(sz, sz, sz, sz, sz, sz);
  REAL*8 acc
  integer arr_len;


  arr_len = sz*sz*sz*sz*sz*sz
#ifdef __INTEL_COMPILER
  write(*, *) "using ifort"
#endif


  ! fill arr1 with random numbers, arr2 with zeros
  call fillrand(arr1, arr_len)
  call fill(arr2, arr_len, 0.0)

  ! do calulcation
  call cpu_time(tstart)
  call outer_func2(arr1, arr2, sz)
  call cpu_time(tfinish)

!  write(*, *) "tstart = ", tstart, ", tfinish = ", tfinish
  print '("outer_func2 time = ", f8.6," second.")', tfinish-tstart


  ! fill arr1 with random numbers, arr2 with zeros
  call fillrand(arr1, arr_len)
  call fill(arr2, arr_len, 0.0)

  ! do calulcation
  call cpu_time(tstart)
  call outer_func2(arr1, arr2, sz)
  call cpu_time(tfinish)

!  write(*, *) "tstart = ", tstart, ", tfinish = ", tfinish
  print '("outer_func2 time = ", f8.6," second.")', tfinish-tstart


STOP
END

! fill array with specified value
SUBROUTINE fill(arr, arr_len, val)
  IMPLICIT NONE
  integer, intent(in) :: arr_len
  REAL*8, intent(inout) :: arr(1:arr_len)
  REAL*8, intent(in) :: val
  integer i
  
  DO i = 1, arr_len
    arr(i) = val
  END DO

END SUBROUTINE

! fill array with random values
SUBROUTINE fillrand(arr, arr_len)
#ifdef __INTEL_COMPILER
  USE IFPORT
#endif

  IMPLICIT NONE
  integer, intent(in) :: arr_len
  REAL*8, intent(inout), DIMENSION(arr_len) :: arr

  ! local variables
  integer i

  DO i = 1, arr_len
    arr(i) = rand()
!    arr(i) = 7
  END DO

END SUBROUTINE

SUBROUTINE outer_func2(u_i, u_ip1, sz)
  IMPLICIT NONE
  integer, intent(in) :: sz
  REAL*8, intent(in), DIMENSION(sz, sz, sz, sz, sz, sz) :: u_i
  REAL*8, intent(inout), DIMENSION(sz, sz, sz, sz, sz, sz) :: u_ip1

  ! local variables
  integer d1, d2, d3, d4, d5, d6, ia, ib, nghost
  REAL*8 delta_12, delta_22, delta_32, delta_42, delta_52, delta_62
  REAL*8 u11, u22, u33, u44, u55, u66

  nghost = 2
  ia = nghost + 1
  ib = sz - nghost

  DO d6=ia, ib
    DO d5=ia, ib
      DO d4=ia, ib
        DO d3 = ia, ib
          DO d2= ia, ib
            DO d1 = ia, ib


              delta_12 = 0.5
              delta_22 = delta_12
              delta_32 = delta_12
              delta_42 = delta_12
              delta_52 = delta_12
              delta_62 = delta_12
             
              u11 = (1.0)*u_i( d1 - 2, d2, d3, d4, d5, d6) + &
                    (2.0)*u_i( d1 - 1, d2, d3, d4, d5, d6) + &
                    (3.0)*u_i( d1    , d2, d3, d4, d5, d6) + &
                    (4.0)*u_i( d1 + 1, d2, d3, d4, d5, d6) + &
                    (5.0)*u_i( d1 + 2, d2, d3, d4, d5, d6)

              u22 = (1.0)*u_i( d1, d2 - 2, d3, d4, d5, d6) + &
                    (2.0)*u_i( d1, d2 - 1, d3, d4, d5, d6) + &
                    (3.0)*u_i( d1, d2    , d3, d4, d5, d6) + &
                    (4.0)*u_i( d1, d2 + 1, d3, d4, d5, d6) + &
                    (5.0)*u_i( d1, d2 + 2, d3, d4, d5, d6)

              u33 = (1.0)*u_i( d1, d2, d3 - 2, d4, d5, d6) + &
                    (2.0)*u_i( d1, d2, d3 - 1, d4, d5, d6) + &
                    (3.0)*u_i( d1, d2, d3    , d4, d5, d6) + &
                    (4.0)*u_i( d1, d2, d3 + 1, d4, d5, d6) + &
                    (5.0)*u_i( d1, d2, d3 + 2, d4, d5, d6)

              u44 = (1.0)*u_i( d1, d2, d3, d4 - 2, d5, d6) + &
                    (2.0)*u_i( d1, d2, d3, d4 - 1, d5, d6) + &
                    (3.0)*u_i( d1, d2, d3, d4    , d5, d6) + &
                    (4.0)*u_i( d1, d2, d3, d4 + 1, d5, d6) + &
                    (5.0)*u_i( d1, d2, d3, d4 + 2, d5, d6)

              u55 = (1.0)*u_i( d1, d2, d3, d4, d5 - 2, d6) + &
                    (2.0)*u_i( d1, d2, d3, d4, d5 - 1, d6) + &
                    (3.0)*u_i( d1, d2, d3, d4, d5    , d6) + &
                    (4.0)*u_i( d1, d2, d3, d4, d5 + 1, d6) + &
                    (5.0)*u_i( d1, d2, d3, d4, d5 + 2, d6)

              u66 = (1.0)*u_i( d1, d2, d3, d4, d5, d6 - 2) + &
                    (2.0)*u_i( d1, d2, d3, d4, d5, d6 - 1) + &
                    (3.0)*u_i( d1, d2, d3, d4, d5, d6    ) + &
                    (4.0)*u_i( d1, d2, d3, d4, d5, d6 + 1) + &
                    (5.0)*u_i( d1, d2, d3, d4, d5, d6 + 2)

              u_ip1( d1, d2, d3, d4, d5, d6) = delta_12*u11 + delta_22*u22 + delta_32*u33 + &
                                               delta_42*u44 + delta_52*u55 + delta_62*u66


            END DO
          END DO
        END DO
      END DO
    END DO
  END DO





END SUBROUTINE

