#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../inc/phi_functs.h"


//PHI Function Struct Defintion.
struct phi_func {
    double *hat;
    double *hat_prime;
};

/////////////////////////////////////////////////////////////////////////////
//Initializing Struct of PHI function arrays to 0.
void InitBasis(struct phi_func *phi, unsigned int size, unsigned int res) {
   unsigned int i = 0;
   unsigned int j = 0;
   for(i = 0; i <= size; i++) { 
      for(j = 0; j <= res; j++) {
         phi[i].hat[j] = 0.00; 
         phi[i].hat_prime[j] = 0.00;
      }
   }
}
////////////////////////////////////////////////////////////////////////////
void SetupBasis(double *x, double *nodes, struct phi_func *phi, unsigned int domain_res, unsigned int num_nodes, unsigned int phi_res) 
{
   unsigned int i = 0;
   unsigned int j = 0;
   double ph = nodes[1] - nodes[0];

   //Filling Right Side of PHI function for phi[0].   
   for(j = 0; j <= domain_res; j++) {
      if((nodes[0] <= x[j]) && (x[j] <= nodes[1])) {
          phi[0].hat[j]         = (nodes[1] - x[j])/ph;
          phi[0].hat_prime[j]   = (-1.00)/ph;
          //fprintf(fp, "%0.3f, %0.3f, %0.3f\n", x[j], phi[0].hat[j], phi[0].hat_prime[j]);
       }
   }
   //Computing PHI Function between the boundary.
   for(i = 1; i <= (num_nodes-1); i++) {
      for(j = 0; j <= domain_res; j++) {
         if ((nodes[i-1] <= x[j]) && (x[j] <= nodes[i])) {
             phi[i].hat[j]       = (x[j] - nodes[i-1])/ph;
             phi[i].hat_prime[j] = 1.00/ph;
             //fprintf(fp, "%0.3f, %0.3f, %0.3f\n", x[j], phi[i].hat[j], phi[i].hat_prime[j]);
         }
         else if ((nodes[i] < x[j]) && (x[j] <= nodes[i+1])) {
             phi[i].hat[j]       = (nodes[i+1] - x[j])/ph;
             phi[i].hat_prime[j] = (-1.00)/ph;
             //fprintf(fp, "%0.3f, %0.3f, %0.3f\n", x[j], phi[i].hat[j], phi[i].hat_prime[j]);
         }
      }
   }
   //Filling Left Side for phi[num_nodes+1] to satisfy boundary condition
   for(j = 0; j <= domain_res; j++) {
      if((nodes[num_nodes-1] <= x[j]) && (x[j] <= nodes[num_nodes])) {
         phi[num_nodes].hat[j]       = ((x[j] - nodes[num_nodes-1])/(ph));
         phi[num_nodes].hat_prime[j] = (1.00)/ph;
         //fprintf(fp, "%0.3f, %0.3f, %0.3f\n", x[j], phi[NUM_NODES].hat[j], phi[NUM_NODES].hat_prime[j]);
      }
   }
}
///////////////////////////////////////////////////////////////////////////

