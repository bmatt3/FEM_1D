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
#define NUM_NODES   500 //Will be Evaluated as N by the code
#define M       (NUM_NODES-1)
#define N       (NUM_NODES-1)

//Code Return Values.
#define FAILURE_LAPACKE     -1
#define FAILURE_BAD_EXIT    -2
#define SUCCESS              0

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

//Preprocessor Matrix Util Function Declarations.
void PrintVector(double *A, unsigned int m);

#ifdef DEBUG
void PrintMatrix(double *A, unsigned int m, unsigned int n);
void PrintVector(double *A, unsigned int m);
#endif 

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
 
  //Loop Counter Variable.
  int i = 0;

  //Sum of all bytes of memory used by variables.
  size_t mem = sizeof(nodes) + sizeof(fx) + sizeof(lv) + sizeof(stiff);
  mem += sizeof(ipiv) + 5*sizeof(lapack_int) + 2*sizeof(double); 
  mem += 2*sizeof(time_t) + sizeof(int);
  
  InitNodes(nodes, a, b, M+1); 
  f(nodes, fx, M+1); 
  StiffMatrix_Assembler(stiff, nodes, M, N);
  Load_Vector_Assembler(fx, nodes, lv, M);

  //Timing Double Precision Symmetric Positive Definite Matrix Solver.
  start = clock();
    info = LAPACKE_dposv(LAPACK_COL_MAJOR, 'U', n, nrhs, stiff, lda, lv, ldb);
  end = clock();
  
  //Checking Return Value. 
  if (info != 0) {
     fprintf(stderr,"LAPACK ERROR : ( %d ) ", info);
     if (info < 0) {
       fprintf(stderr, " illegal value of one or more");
       fprintf(stderr, "arguments -- no computation performed.\n");
     }
     else {
       fprintf(stderr, " failure in the course of computation.\n");
     }
     return FAILURE_LAPACKE;
  } 
  else {
     //On success we write to the solution to disk for post processing.
     FILE *fp = fopen("./data/Solution.dat", "w");
     fprintf(fp, "%f, %f\n",  0.00, 0.00);
     for(i = 0; i < M; i++) {
       fprintf(fp, "%f, %f\n",  nodes[i+1], lv[i]);
     }   
     fprintf(fp, "%f, %f\n",  1.00, 0.00);
     fclose(fp); 

     //Printing Total Memory Usage By Program
     if (mem < pow(1024,2)) {
        printf("Program Memory Utilization = %0.3f (KiB)\n", (double) (mem/(1024)));
     }
     else if (mem < pow(1024,3)) {
        printf("Program Memory Utilization = %0.3f (MiB)\n", (double) (mem/(pow(1024,2))));
     }
     else {
        printf("Program Memory Utilization = %0.3f (GiB)\n", (double) (mem/pow(1024,3)));
     }

     //Printing Symmetric Positive Matrix Algorithm Runtime.
     printf("dposv_run time = %f\n\n", ((double) (end - start)) / CLOCKS_PER_SEC);
     return SUCCESS;
  }
  //We should never get here. But we will error handle for safety.
  return FAILURE_BAD_EXIT;
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

    #ifdef DEBUG 
    printf("Nodal Values:");
    PrintVector(nodes, num_nodes);
    printf("\n\n");
    #endif 
}
/////////////////////////////////////////////////////////////////////////////
//Using the Galerkin Method to assemble the Mass Matrix.
//Using Simpson Method for Integration Scheme.
void StiffMatrix_Assembler(double *stiff, double *nodes, unsigned int rows, unsigned int cols)
{
  unsigned int i = 0;
  double       h = nodes[1] - nodes[0];
  
  FillMatrixDiagonal(stiff, rows, cols,  0, (2.0/h));
  FillMatrixDiagonal(stiff, rows, cols,  1, (-1.0/h));
  FillMatrixDiagonal(stiff, rows, cols, -1, (-1.0/h));

  #ifdef DEBUG
  printf("Printing Stiffness Matrix Values (Column Major Order):\n");
  PrintMatrix(stiff, rows, cols);
  printf("\n\n");
  #endif
}
///////////////////////////////////////////////////////////////////////////
//Using the Galerkin Method to Assemble to Load Vector
void Load_Vector_Assembler(double *fx, double *nodes, double *lv, \
     unsigned int num_nodes)
{
  unsigned int i = 1;
  double       h = nodes[1] - nodes[0];

  MATRIX_ELEMENT(lv,num_nodes,1,0,0) = MATRIX_ELEMENT(lv,num_nodes,1,0,0) + (fx[1]*(h/2));
  for(i = 1; i <= (num_nodes-2); i++) {
     MATRIX_ELEMENT(lv,num_nodes,1,i,0) = MATRIX_ELEMENT(lv,num_nodes,1,i,0) + (fx[i+1]*(h));
  }
  MATRIX_ELEMENT(lv, num_nodes, 1, num_nodes-1,0) =  MATRIX_ELEMENT(lv,num_nodes,1,num_nodes-1,0) + fx[num_nodes]*(h/2);

  #ifdef DEBUG
  printf("Load Vector Values:\n");
  PrintVector(lv, num_nodes);
  printf("\n\n");
  #endif
}
//////////////////////////////////////////////////////////////////////////
//Evaluating f(x) at the nodal points.
double f(double *nodes, double *fx, unsigned int num_nodes) {
   unsigned int i = 0;
   for(i = 0; i <= (num_nodes-1); i++) {
       fx[i] = nodes[i]*sin(nodes[i]);
   }

   #ifdef DEBUG
   printf("f(x) values at nodal points: ");
   PrintVector(fx, num_nodes);
   printf("\n\n");
   #endif 
}
