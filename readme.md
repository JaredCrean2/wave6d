This repository contains finite difference solvers for the wave equation
in 1 through 6 dimensions.  An explicit Runge-Kutta method is used for
the time discretization and finite differences are used for the spatial
discretization.  A script is used to generate the spatial kernel, and can
use an arbitrary order stencil in each direction.

## Code Structure
### Code Generation
The `generate_code.jl` script will generate several source files and put them
in the `src/generated` directory (make sure the directory exists first).
There are several parameters at the top of `generate_code.jl` that control
the maximum dimension solver to generate and the stencil to use.  If these are
modified the corresponding includes at the top of `WaveSolver.jl` must be
modified as well.  The naming scheme of the files is fairly simple.

### The Solver
`WaveSolver.jl` defines the `WaveSolver` module, which contains the solver.
The only exported method is `runcase(fname)`, where `fname` is the name of 
an input file.  The format of the input file is documented in the function
`parseinput`.  The file `generate_input.jl` contains a function that assists
in the creation of input files.

The solver has 3 modes of operation: no loop blocking, partial loop blocking,
and Hilbert curve.  In no loop blocking mode, the function that evaluates the
spatial derivatives loops over the array linearly.  In partial loop blocking
mode, a number of loops (specified by the input file) are blocked, where the
block size is specified by the input file.  Note that the number of points
in the grid must evenly divide into the block size (the solver currently 
declines to implement the cleanup loops required to handle the non-evenly
divisible case).  In Hilbert curve mode, the points on the grid are traversed
according to a Hilbert curve.  The array dimensions must all be the same and
be powers of 2 in this mode.

The `runconvergence.jl` script can be used to run converence studies.
Several parameters at the top of the file can be used to set the conditions
of the run.


### Parallelization
The code is parallelized using MPI.  The parallel initialization code 
assigns the MPI processes to a grid that is as close to square as possible,
using an efficent factorization based algorithm.  Although I cannot formally
prove optimality, I have strong reason to believe the algorithm is optimal and
it agrees with an exhaustive search in all cases tested.

One the MPI grid is defined, points are assigned to each processor. 
The points are assigned using a median-remainder scheme where some number
of points are assigned to the first N processors along an axis, and the
remainder are assigned to the last process.  Given this constraint,
this point distributions it produces are fairly close to
optimal (there are a few corner cases where sub-optimal decompositions are
produced in order to avoid degenerate cases.  These mostly occur grids
that are very skewed in one dimension, or have very small numbers of points).
Also note that the MPI grid is defined indendently of the number of points
in each direction.  This enables decoupling the decomposition of the MPI ranks
and the decomposition of the node points (although I think the algorithm can be
adapted to couple the decompositions without too much difficulty).

Note that it is the callers responsibility to initialize MPI before calling
`runcase`, and to finalize it afterwards

Also, there is a strange behavior with MPI on some platforms where `MPI.Wait`
does not behave correctly when operation in a previous send operation.  
MPI returns a nonsense error value and the waiting does not occur.  Waiting
on a receive operation appears to work correctly.  There is a global const
flag named `FORCE_SYNC`.  This inserts a `MPI.Barrier` into each function
evaluation which prevents sending data to another process before it received
the previous data.

### Time Stepping
Either the classical 4th order explict Runge-Kutta method or the 5 stage, 4th
order low storage varient by Carpenter and Kennedy can be used.  Tests show
the low storage varient to be 25% slower than the classical version for the
same CFL, exactly as expected.  The flag `USE_LOW_STORAGE` at the top of
`WaveSolver.jl` controls which method is used.





## Performance Testing

### Comparison with Fortran

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

### Systems of Equations
Solving a `d` dimensional system of equations requires an array with `d+1`
dimenions. The question is whether the dimension corresponding to the equation
number should be the first or the last.  The file `test_inline6.jl` contains
a test for a 5 dimensional system.  `outer_func` puts the equation index first,
`outer_func2` puts it last.  Results:

```
testing outer_func
  0.361604 seconds
  0.362164 seconds
testing outer_func2
  0.215627 seconds
  0.214228 seconds
```

A look at the llvm code shows `outer_func` scalarized the loop. For absolute
performance, it looks like vectorized instructions won out over spatial
locality.


