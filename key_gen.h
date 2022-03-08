
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>


#include "matrix.h"

#include "param.h"
#include "util.h"


int disjoint_test(gf * u, gf * v );
int Test_disjoint(gf * L,int n);
gf* Random_Vect(int m);
void generate_random_vector(int m, gf *vect);
void init_random_element(gf *U);

gf* Init_Random_U();
void binary_quasi_dyadic_sig(int m, int n, int t, int * b, gf * h_sig, gf * w );
void cauchy_support(gf * Support, gf * W,gf * w);
int key_pair(unsigned char * pk, unsigned char * sk);
void set_y_from_uvz(gf_t * u, gf_t* v, gf_t * Z, gf_t* y);

int Try(binmat_t A);
