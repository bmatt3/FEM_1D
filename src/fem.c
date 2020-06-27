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
#include <stdint.h>
#include <math.h>
#include <lapacke.h>
#include "../inc/matrix_utils.h"

//Defining Constants for Basis Function Construction.
#define NUM_NODES   8
#define DOMAIN_RES  100

//Code Return Values.
#define FAILURE     0
#define SUCCESS     1

//Column Major Order Preprocessor Macros.
//**Not used because we are no longer using dgesv routine.
//**Since Fortran LAPACK is written in Fortran it expects column major order
//**oriented arrays. 

//#define MATRIX_IDX(n,i,j) j*n + i
//#define MATRIX_ELEMENT(A,m,n,i,j) A[MATRIX_IDX(m,i,j)]

///////////////////////////////////////////////////////////////////////////////

//Function Declarations:
void InitDomain(double *x, double a, double b, double n);

void InitNodes(double *nodes, double *x, unsigned int num_nodes, \
     unsigned int x_size);

void StiffMatrix_Assembler(double *nodes, double *dl, \
     double *d, double *du, unsigned int num_nodes);

double f(double *nodes, double *fx, unsigned int num_nodes);

void Load_Vector_Assembler(double *fx, double *nodes, \
     double *lv,unsigned int num_nodes);

//Matrix Util Function Declarations.
void PrintVector(double *A, unsigned int m);

////////////////////////////////////////////////////////////////////////////////
//CODE START:
int main(int argc, char *argv[]) 
{ 

  double  x[DOMAIN_RES+1]   = {0.00}; //X Domain Values.
  double fx[NUM_NODES-1]    = {0.00}; //f(x) values at nodal points.
  double dl[NUM_NODES-3]    = {0.00}; //Sub diagonal values of Mass Mat.
  double  d[NUM_NODES-2]    = {0.00}; //Diagonal values of Mass Mat.
  double du[NUM_NODES-3]    = {0.00}; //Super diagonal values of Mass Mat.
  double lv[NUM_NODES-1]    = {0.00}; //LoadVector Values.
  double nodes[NUM_NODES-1] = {0.00};

  double a = 0.00;    //X-Domain Start
  double b = 1.00;    //X-Domain End

  //LAPACKE_DGTSV argument Initialization.
  int n    = NUM_NODES-1; 
  int lda  = NUM_NODES-1; 
  int nrhs = 1;
  int ldb  = 1;
  int info = 0;

  
  InitDomain(x, a, b, DOMAIN_RES);
  InitNodes(nodes, x, NUM_NODES-1, DOMAIN_RES); 
  f(nodes, fx, NUM_NODES-1);
  StiffMatrix_Assembler(nodes, dl, d, du, NUM_NODES-1);
  Load_Vector_Assembler(fx, nodes, lv, NUM_NODES-1);
  
  //printf("X Vector:\n");
  //PrintVector(x, DOMAIN_RES+1);
  //printf("\n\n");
 
  printf("Node Vector:\n");
  PrintVector(nodes, NUM_NODES-2);
  printf("\n\n"); 

  printf("Fx Vector:\n");
  PrintVector(fx, NUM_NODES-1);
  printf("\n\n");
/*
  printf("Upper Diagonal Vector:\n");
  PrintVector(du, NUM_NODES-2);
  printf("\n\n");
  
  printf("Diagonal Vector:\n");
  PrintVector(d, NUM_NODES-1);
  printf("\n\n");

  printf("Lower Diagonal Vector:\n");
  PrintVector(dl, NUM_NODES-2);
  printf("\n\n");
  
  printf("LoadVector:\n");
  PrintVector(lv, NUM_NODES-1);
  printf("\n\n");
  */

  //Calling LAPACKE_dgtsv solves triangular matrix in O(n) time.
  info = LAPACKE_dgtsv(LAPACK_ROW_MAJOR, n, nrhs, dl, d, dl, lv, ldb);

  if (info != 0) {
     fprintf(stderr,"LAPACK ERROR : ( %d ) ", info);
     if (info < 0) {
       fprintf(stderr, " illegal value of one or more");
       fprintf(stderr, "arguments -- no computation performed.\n");
     }
     else {
       fprintf(stderr, " failure in the course of computation.\n");
     }
     return FAILURE;
  } 
  else {
     printf("\n\nSolution:\n");
     PrintVector(lv, NUM_NODES-2);
     //FILE *fp = stdout; 
     FILE *fp = fopen("./data/Solution.dat", "w");
     int i = 0;
     //fprintf(fp, "%f, %f\n",  0.00, 0.00);
     for(i = 0; i < (NUM_NODES-2); i++) {
       fprintf(fp, "%f, %f\n",  nodes[i], lv[i]);
        printf( "%f, %f\n",  nodes[i], lv[i]);
     }   
     //fprintf(fp, "%f, %f\n",  1.00, 0.00);
     fclose(fp);
  }
  return SUCCESS;
}
//////////////////////////////////////////////////////////////////////////////
//Initialzing the problem domain.
void InitDomain(double *x, double a, double b, double n) 
{
    unsigned int i   = 1;
    double       h   = (b - a)/n;

    x[0] = a;
    for(i=1; i <= (DOMAIN_RES-1); i++) {
       x[i] = x[i-1] + h;
    }
    x[DOMAIN_RES] = b;
}
//////////////////////////////////////////////////////////////////////////////
//Computing uniformly spaced nodal values on the x-axis,
//And allocating them to the nodes array.
void InitNodes(double *nodes, double *x, unsigned int num_nodes, \
     unsigned int x_size) 
{
    unsigned int i = 0;   
    double h = (1.00 - 0.00)/(num_nodes);
    for(i = 1; i <= num_nodes; i++) {
       nodes[i] = nodes[i-1] + h;
    }
}
/////////////////////////////////////////////////////////////////////////////
//Using the Galerkin Method to assemble the Mass Matrix.
//Using Simpson Method for Integration Scheme.
void StiffMatrix_Assembler(double *nodes, double *dl, \
     double *d, double *du, unsigned int num_nodes)
{
  unsigned int i = 0;
  double       h = nodes[2] - nodes[1];  
  
  //Filling Lower, Upper, and Central Diagonal 
  for(i = 0; i <= (num_nodes-2); i++) {
     dl[i] = (-1.0)/h;
     du[i] = (-1.0)/h;  
     d[i]  = (2.0)/h;
  }
  d[num_nodes-1] = (2.0)/h;
}
///////////////////////////////////////////////////////////////////////////
//Using the Galerkin Method to Assemble to Load Vector
void Load_Vector_Assembler(double *fx, double *nodes, double *lv, \
     unsigned int num_nodes)
{
    unsigned int i = 1;
    double       h = nodes[2] - nodes[1];

    lv[0] = lv[0] + (fx[1]*(h/2));
    for(i = 1; i <= (num_nodes-3); i++) { 
       lv[i] = lv[i] + (fx[i+1]*(h));
    }
    lv[num_nodes-2] = lv[num_nodes-2] + (fx[num_nodes-1]*(h/2));
}
//////////////////////////////////////////////////////////////////////////
//Evaluating f(x) at the nodal points.
double f(double *nodes, double *fx, \
       unsigned int num_nodes) 
{
   unsigned int i = 0;
   fx[0] = 1.00;
   for(i = 1; i <= num_nodes; i++) {
       fx[i] = 1.00; 
   }
}
