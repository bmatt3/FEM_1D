#ifdef  PHI_FUNCT_GUARD
#define PHI_FUNCT_GUARD

struct phi_func {
    double *hat;
    double *hat_prime;
};

void InitBasis(struct phi_func *phi, unsigned int size, unsigned int res);

void SetupBasis(double *x, double *nodes, struct phi_func *phi, unsigned int domain_res, \
		unsigned int num_nodes, unsigned int phi_res) ;

#endif
