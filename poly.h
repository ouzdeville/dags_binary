
#ifndef POLY_H
#define POLY_H

#include "gf.h"
#include "param.h"

typedef struct polynome {
  int deg, size;
  gf_t * coeff;
} * poly_t; /* polynomial has coefficients in the finite field */

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

#define poly_deg(p) ((p)->deg)
#define poly_size(p) ((p)->size)
#define poly_set_deg(p, d) ((p)->deg = (d))
#define poly_coeff(p, i) ((p)->coeff[i])
#define poly_set_coeff(p, i, a) ((p)->coeff[i] = (a))
#define poly_addto_coeff(p, i, a) ((p)->coeff[i] = gf_add((p)->coeff[i], (a)))
#define poly_multo_coeff(p, i, a) ((p)->coeff[i] = gf_mul((p)->coeff[i], (a)))
#define poly_tete(p) ((p)->coeff[(p)->deg])

/****** poly.c ******/

int poly_calcule_deg(poly_t p);
poly_t poly_alloc(int d);
poly_t poly_alloc_from_string(int d, const unsigned char * s);
poly_t poly_copy(poly_t p);
void poly_free(poly_t p);
void poly_set_to_zero(poly_t p);
void poly_set(poly_t p, poly_t q);
poly_t poly_mul(poly_t p, poly_t q);
void poly_rem(poly_t p, poly_t g);
void poly_sqmod_init(poly_t g, poly_t * sq);
void poly_sqmod(poly_t res, poly_t p, poly_t * sq, int d);
poly_t poly_gcd(poly_t p1, poly_t p2);
poly_t poly_quo(poly_t p, poly_t d);
gf_t poly_eval(poly_t p, gf_t a);
int poly_degppf(poly_t g);
void poly_eeaux(poly_t * u, poly_t * v, poly_t p, poly_t g, int t);

poly_t * poly_syndrome_init(poly_t generator, gf_t *support, int n);
poly_t * poly_sqrtmod_init(poly_t g);
poly_t poly_randgen_irred(int t, int (*u8rnd)());
void poly_add_free (poly_t r, poly_t a, poly_t b);

void Sq_mod(poly_t p, poly_t * sq, int t );
void poly_display(poly_t p, int n);
void Tr_pol(poly_t * tr,poly_t p, int e, int m) ;
int BTA_round(poly_t sigma, gf_t * res, poly_t * tr, int e, int d, int num_tr) ;
int BTA(poly_t sigma, gf_t * res) ;
#endif /* POLY_H */
