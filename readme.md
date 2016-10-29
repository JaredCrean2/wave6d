This repository contains finite difference solvers for the wave equation
in 1 through 6 dimensions.  An explicit Runge-Kutta method is used for
the time discretization and finite differences are used for the spatial
discretization.  A script is used to generate the spatial kernel, and can
use an arbitrary order stencil in each direction.

## Performance Testing

`src/test_inline.jl` and `src/test_inline.f90` have test implementations of a 
6 dimensional kernel.  The performance results are

Fortran:
```
  outer_func2 time = 0.160010 second.
```
Julia:
```
testing outer_func
  0.159087 seconds
  0.157757 seconds
testing outer_func2
  0.150705 seconds
  0.151781 seconds
```

The disctinction between `outer_func` and `outer_func2` is that `outer_func`
passes a tuple of indices to the kernel function, which unpacks them and does
the computation for a single grid point.  `outer_func2` contains both the
loops and the code for the kernal in a single function.  Only `outer_func2`
was implemented in Fortran.  A look at the assembly code shows both languages
vectorized the computation.  The conclusions from this data are that Julia
is as fast as Fortran for this kind of computation and that separating
the kernel into a separate function still gives good performance.

