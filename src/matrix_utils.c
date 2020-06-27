#include <stdio.h>
#include <stdlib.h>
#include "../inc/matrix_utils.h"

#define MATRIX_IDX(n,i,j) j*n + i
#define MATRIX_ELEMENT(A,m,n,i,j) A[MATRIX_IDX(m,i,j)]

void PrintMatrix(double *A, unsigned int m, unsigned int n) {
  unsigned int i = 0;
  unsigned int j = 0;
  for(i = 0; i < m; i++) {
     for(j = 0; j < n; j++) {
        printf("%8.4f", MATRIX_ELEMENT(A, m, n, i, j));
      }
      printf("\n");
   }
}

void FillMatrixDiagonal(double *A, unsigned int m, unsigned int n, unsigned int k, double val) {
  unsigned int i = 0;
  unsigned int j = 0;
  for(i = 0; i < n; i++) {
    for(j = 0; j < m; j++) {
       if(i == (j+k)) {
          MATRIX_ELEMENT(A,m,n,i,j) = val;
       }
    }
  }
}

void PrintVector(double *A, unsigned int m) {
   unsigned int i = 0;
   for(i = 1; i <= m; i++) {
      printf("%d \t %f\n", i-1, A[i-1]);
   }
   printf("\n");
}


