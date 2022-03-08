/******************************************************************************************
* BINARY DAGS: Key Encapsulation using binary Dyadic GS Codes.                            *
* This is the binary version of DAGS  submitted to the NIST Post=Quantum Cryptography.    *
* For any other usage , contact the author(s) to ascode_dimmension permission.            *
*******************************************************************************************
*******************************************************************************************
*********************** MACROs for certain operations ************************************/



#ifndef GF_H
#define GF_H

typedef unsigned short gf_t;
typedef unsigned short gf;

int gf_extension_degree, gf_cardinality, gf_multiplicative_order;
gf_t * gf_log;
gf_t * gf_exp;





#define gf_extd() gf_extension_degree
#define gf_card() gf_cardinality
#define gf_ord() gf_multiplicative_order

#define gf_unit() 1
#define gf_zero() 0
#define gf_add(x, y) ((x) ^ (y))
#define gf_exp(i) gf_exp[i] /* alpha^i */
#define gf_log(x) gf_log[x] /* return i when x=alpha^i */

// residual modulo q-1
// when -q < d < 0, we get (q-1+d)
// when 0 <= d < q, we get (d)
// when q <= d < 2q-1, we get (d-q+1)
#define _gf_modq_1(d) (((d) & gf_ord()) + ((d) >> gf_extd()))
/* we obtain a value between 0 and (q-1) included, the class of 0 is
represented by 0 or q-1 (this is why we write _K->exp[q-1]=_K->exp[0]=1)*/

#define gf_mul_fast(x, y) ((y) ? gf_exp[_gf_modq_1(gf_log[x] + gf_log[y])] : 0)
#define gf_mul(x, y) ((x) ? gf_mul_fast(x, y) : 0)
#define gf_square(x) ((x) ? gf_exp[_gf_modq_1(gf_log[x] << 1)] : 0)
#define gf_sqrt(x) ((x) ? gf_exp[_gf_modq_1(gf_log[x] << (gf_extd()-1))] : 0)
// Try to devide by zero and get what you deserve!
#define gf_div(x, y) ((x) ? gf_exp[_gf_modq_1(gf_log[x] - gf_log[y])] : 0)
#define gf_inv(x) gf_exp[gf_ord() - gf_log[x]]

/****** gf.c ******/

int gf_init(int extdeg);
gf_t gf_rand(int (*u8rnd)());
gf_t gf_pow(gf_t x, int i);

#endif /* GF_H */
