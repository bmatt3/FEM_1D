////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//Problem Description:                                                        //
//The problem that we are solving is -u''(x) = f(x).                          //
//Reference --> "An Introduction to the Finite Element Method                 //
//For Differential Equations" (Mohammad Asadzadeh, Jan 20, 2010).             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>
#include <time.h>
#include "../inc/matrix_utils.h"

//Defining Constants for Basis Function Construction.
#define NUM_NODES   7 //Will be Evaluated as N by the code
#define M       (NUM_NODES-1)
#define N       (NUM_NODES-1)

//Code Return Values.
#define FAILURE     1
#define SUCCESS     0

//Defining Column Major Order Preprocessor Macros.
#define MATRIX_IDX(n,i,j) j*n + i
#define MATRIX_ELEMENT(A,m,n,i,j) A[MATRIX_IDX(m,i,j)]

///////////////////////////////////////////////////////////////////////////////
//Preprocessor Function Declarations:
void InitNodes(double *nodes, double a, double b, unsigned int num_nodes);

double f(double *nodes, double *fx, unsigned int num_nodes);

void StiffMatrix_Assembler(double *stiff, double *nodes, unsigned int rows, unsigned int cols);

void Load_Vector_Assembler(double *fx, double *nodes, \
     double *lv,unsigned int num_nodes);

//Preprocessro Matrix Util Function Declarations.
void PrintVector(double *A, unsigned int m);

//void PrintMatrix(double *A, unsigned int m, unsigned int n);

void FillMatrixDiagonal(double *A, unsigned int m, unsigned int n, \
     unsigned int k, double val);

////////////////////////////////////////////////////////////////////////////////

//CODE START:
int main(int argc, char *argv[]) 
{ 
  double a           = 0.00; //X-Domain Start
  double b           = 1.00; //X-Domain End
  double nodes[M+1]  = {0};  //Nodal points on the x-axis.
  double fx[M+1]     = {0};  //f(x) values at nodal points.
  double lv[M]       = {0};  //LoadVector Values.
  double stiff[M*N]  = {0};  //Stiffness Matrix.

  lapack_int n       = M;   //Order of Matrix A.
  lapack_int lda     = M;   //Leading Dimension Matrix A.
  lapack_int ipiv[M] = {0}; //Pivot Array. Entries Initialized to 0.
  lapack_int nrhs    = 1;   //Number of right hand sides.
  lapack_int ldb     = M;   //Leading Dimension of B.
  lapack_int info    = 0;   //Return Value.
  
  //Timing Variables 
  time_t start; 
  time_t end; 
 
  size_t memory_usage = sizeof(nodes) + sizeof(fx) + sizeof(lv) + sizeof(stiff);
  memory_usage += sizeof(ipiv) + 6*sizeof(int) + 2*sizeof(double) + sizeof(size_t);
  memory_usage += 2*sizeof(clock_t);

  printf("RAM USAGE = %0.1f (KB)\n",  memory_usage/1.0e3);
  printf("RAM USAGE = %0.1f (MB)\n", memory_usage/1.0e6);
  printf("RAM USAGE = %0.1f (GB)\n", memory_usage/1.0e9);
  
  InitNodes(nodes, a, b, M+1); 
  f(nodes, fx, M+1); 
  StiffMatrix_Assembler(stiff, nodes, M, N);
  Load_Vector_Assembler(fx, nodes, lv, M);

  //Double Precision, Generalized Matrix Solver
  //start = clock();
  //  info = LAPACKE_dgesv(LAPACK_COL_MAJOR, n, nrhs, stiff, lda, ipiv, lv, ldb);
  //end = clock();
  //printf("dgesv_run time = %f\n\n", ((double) (end - start)) / CLOCKS_PER_SEC);
  
  //Double Precision, Symmetric Positive Definite Matrix Solver.  
  //Cholskey Decomposition + Linear Solve =( (1/3)*O(n^3) + O(n^2)) Runtime)
  start = clock();
    info = LAPACKE_dposv(LAPACK_COL_MAJOR, 'U', n, nrhs, stiff, lda, lv, ldb);
  end = clock();
  printf("dposv_run time = %f\n\n", ((double) (end - start)) / CLOCKS_PER_SEC);

  if (info != 0) {
     fprintf(stderr,"LAPACK ERROR : ( %d ) ", info);
     if (info < 0) {
       fprintf(stderr, " illegal value of one or more");
       fprintf(stderr, "arguments -- no computation performed.\n");
     }
     else {
       fprintf(stderr, " failure in the course of computation.\n");
     }
  } 
  else {
     FILE *fp = fopen("./data/Solution.dat", "w");
     fprintf(fp, "%f, %f\n",  0.00, 0.00);
     for(info = 0; info < M; info++) {
       fprintf(fp, "%f, %f\n",  nodes[info+1], lv[info]);
     }   
     fprintf(fp, "%f, %f\n",  1.00, 0.00);
     fclose(fp); 
  }
  return SUCCESS;
}
//////////////////////////////////////////////////////////////////////////////
//Computing uniformly spaced nodal values on the x-axis,
//And allocating them to the nodes array.
void InitNodes(double *nodes, double a, double b, unsigned int num_nodes) {
    unsigned int i = 0;   
    double h = (b-a)/((double) num_nodes);
    nodes[0] = a;
    for(i = 1; i <= num_nodes+1; i++) {
       nodes[i] = nodes[i-1] + h;
    }
}
/////////////////////////////////////////////////////////////////////////////
//Using the Galerkin Method to assemble the Mass Matrix.
//Using Simpson Method for Integration Scheme.
void StiffMatrix_Assembler(double *stiff, double *nodes, unsigned int rows, unsigned int cols)
{
  unsigned int i = 0;
  double       h = nodes[1] - nodes[0];
  FillMatrixDiagonal(stiff, rows, cols,  0, (2.0/h));   //Writing To Main Diagonal.
  FillMatrixDiagonal(stiff, rows, cols,  1, (-1.0/h));  //Writing To 1st Super-Diagonal.
  FillMatrixDiagonal(stiff, rows, cols, -1, (-1.0/h));  //Writing To 1st Sub-Diagonal.
}
///////////////////////////////////////////////////////////////////////////
//Using the Galerkin Method to Assemble to Load Vector
void Load_Vector_Assembler(double *fx, double *nodes, double *lv, \
     unsigned int num_nodes)
{
    unsigned int i = 1;
    double       h = nodes[1] - nodes[0];
     
    //Filling Load Vector in Column Major Order.
    MATRIX_ELEMENT(lv,num_nodes,1,0,0) = MATRIX_ELEMENT(lv,num_nodes,1,0,0) + (fx[1]*(h/2));
    for(i = 1; i <= (num_nodes-2); i++) {
       MATRIX_ELEMENT(lv,num_nodes,1,i,0) = MATRIX_ELEMENT(lv,num_nodes,1,i,0) + (fx[i+1]*(h));
    }
    MATRIX_ELEMENT(lv, num_nodes, 1, num_nodes-1,0) = MATRIX_ELEMENT(lv,num_nodes,1,num_nodes-1,0) + fx[num_nodes]*(h/2);
}
//////////////////////////////////////////////////////////////////////////
//Evaluating f(x) at the nodal points.
double f(double *nodes, double *fx, unsigned int num_nodes) {
   unsigned int i = 0;
   for(i = 0; i <= (num_nodes-1); i++) {
       fx[i] = 1.00;
   }
}
