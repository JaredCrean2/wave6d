#ifndef __TEST_H__
#define __TEST_H__

#include "hilbert.h"
#include "stdio.h"
#include "stdlib.h"
#include "stdint.h"


void printArray(FILE* f, bitmask_t arr[], unsigned int len);
void printArray16(FILE* f, uint16_t arr[], int len);
int idx(const int i, const int j, const int si, const int sj);

void loadNpoints(int n, int dims, bitmask_t coords[], uint16_t idxs[]);
int checkDimensions(int dims, int npoints);
#endif 
