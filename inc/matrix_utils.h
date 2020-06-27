#ifdef  MATRIX_UTIL_GUARD
#define MATRIX_UTIL_GUARD

#define MATRIX_IDX(n,i,j) j*n + i
#define MATRIX_ELEMENT(A,m,n,i,j) A[MATRIX_IDX(m,i,j)]

void PrintMatrix(double *A, unsigned int m, unsigned int n);
void FillMatrixDiagonal(double *A, unsigned int m, unsigned int n, unsigned int k, double val);
void PrintVector(double *A, unsigned int m);

#endif 
