#include "hilbert.h"
#include "stdio.h"
#include "stdlib.h"
#include "stdint.h"

void printArray(FILE* f, bitmask_t arr[], unsigned int len)
{
  for (unsigned int i=0; i < len; ++i)
    fprintf(f, "%llu ", arr[i]);

  fprintf(f, "\n");
}

void printArray16(FILE* f, uint16_t arr[], int len)
{
  for (unsigned int i=0; i < len; ++i)
    fprintf(f, "%d ", arr[i]);

  fprintf(f, "\n");
}
// column major indexing, 0 based
int idx(const int i, const int j, const int si, const int sj)
{
  return i + j*si;
}

// get the next n points on the curve
// use unit16_t for the indices to save space
// idxs is a dim x n array to be populated with the cartesian indices
// on entry, coords_prev should contain the first point on the current 
// section of the curve.  On exit it will contain the first point on the
// next section of the curve.  Thus the first call should have coords 
// as all zeros and then use teh same coords array for subsequent calls
void loadNpoints(int n, int dims, bitmask_t coords[], uint16_t idxs[])
{
  // integer division
  unsigned int nbits_per_dim = (8*sizeof(bitmask_t) / dims);


  for (int i = 0; i < n; ++i)
  {    // copy to idxs
    for (int j = 0; j < dims; ++j)
      idxs[idx(j, i, dims, n)] = coords[j] + 1;  // convert to 1 based indexing

    hilbert_incr(dims, nbits_per_dim, coords);
  }

}  // end function loadNpoints

// check that there are enough bits per dimension to represent the
// maximum number of points in any dimension
int checkDimensions(int dims, int npoints  )
{
  
  int ret_status = 0;
  unsigned int nbits_per_dim = (8*sizeof(bitmask_t) / dims);
  unsigned long long int max_size = 2;
  for (int i = 1; i < nbits_per_dim; ++i)
    max_size *= 2;

  printf("nbits_per_dim = %d\n", nbits_per_dim);
  printf("max_size = %d\n", max_size);
  printf("npoints = %d\n", npoints);
  if ( max_size < npoints) // -1 becaue zero is representable
  {
    fprintf(stderr, "Error: using %d bits per dimension can represent a maximum of %d points, %d is too many\n", nbits_per_dim, max_size, npoints);
    ret_status = 1;
  }

  return ret_status;
}



int main(int argc, char* argv[])
{

  printf("Hello world\n");

  printf("size of bitmask_t = %d\n", sizeof(bitmask_t));
  printf("size of halfmask_t = %d\n", sizeof(halfmask_t));

  unsigned int ndims = 2; // test 2D hilbert curve
  unsigned int nbits_per_dim = (8*sizeof(bitmask_t) / ndims) - 2;

  bitmask_t coords[ndims];

  // zero out coords
  for (int i=0; i < ndims; ++ i)
    coords[i] = 0;

  FILE* f = fopen("hdata.dat", "w");
//  printf("initially, coords = ");
  printArray(f, coords, ndims);
//  fprintf("\n");


  for (int i = 0; i < 1023; ++i)
  {
    hilbert_incr(ndims, nbits_per_dim, coords);

//    printf("after increment %d, coords = ", i);
    printArray(f, coords, ndims);
  }

  fclose(f);

  int n = 64;
  uint16_t* idxs = (uint16_t *)malloc(2*n*ndims);
  // zero out coords
  for (int i=0; i < ndims; ++i)
    coords[i] = 0;

  loadNpoints(n, ndims, coords, idxs);

  f = fopen("hdata2.dat", "w");
  for (int i=0; i < n; ++i)
    printArray16(f, idxs + idx(0, i, ndims, n), ndims);

  fclose(f);

  int ret_stat = checkDimensions(2, 1024);
  printf("ret_stat = %d\n", ret_stat);

  ret_stat = checkDimensions(4, 65537);
  printf("ret_stat = %d\n", ret_stat);


  return 0;
}

