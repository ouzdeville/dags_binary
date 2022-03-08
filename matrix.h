/******************************************************************************************
 * BINARY DAGS: Key Encapsulation using binary Dyadic GS Codes.                            *
 * This is the binary version of DAGS  submitted to the NIST Post=Quantum Cryptography.    *
 * For any other usage , contact the author(s) to ascode_dimmension permission.            *
 *******************************************************************************************
 */

#include "rng.h"
#include "gf.h"
#include "param.h"
#ifndef _MATRIX_H
#define _MATRIX_H

typedef struct matrix {
	int rown; //number of rows.
	int coln; //number of columns.
	int rwdcnt; //number of words in a row
	int alloc_size; //number of allocated bytes
	unsigned long **elem; //row index.
} binmat_t;

typedef struct Matrix {
	unsigned int rown;
	unsigned int coln;
	int rwdcnt; //number of words in a row
	gf_t **coeff;
} mat_t;

typedef struct _dyadic_Matrix {
	//unsigned int ordre;
	int rwdcnt; //number of words in a row
	int alloc_size; //number of allocated bytes
	unsigned long *coeff;
} D_mat;

typedef struct quasi_dyadic_Matrix {
	unsigned int rown;
	unsigned int coln;
	unsigned int ordre;
	D_mat **dblock;
} QD_mat;

#define mat_coeff(A, i, j) (((A).elem[(i)][(j) / BITS_PER_LONG] >> (j % BITS_PER_LONG)) & 1)
//#define mat_row(A, i) ((A)->elem + ((i) * A->rwdcnt))
#define mat_set_coeff_to_one(A, i, j) ((A).elem[(i)][(j) / BITS_PER_LONG] |= (1UL << ((j) % BITS_PER_LONG)))
#define mat_change_coeff(A, i, j) ((A).elem[(i)][ (j) / BITS_PER_LONG] ^= (1UL << ((j) % BITS_PER_LONG)))
#define swaprows(A, i, k) {unsigned long * pt = A.elem[i]; A.elem[i] =A.elem[k]; A.elem[k] = pt;}
//#define mat_set_to_zero(R) memset((R)->elem,0,(R)->alloc_size);
//#define mat_set_to_zero(R) memset((R)->elem,0,(R)->alloc_size);

//For Dyadic Blocks
#define D_mat_coeff(A, j) (((A).coeff[(j) / BITS_PER_LONG] >> (j % BITS_PER_LONG)) & 1)
#define D_mat_set_coeff_to_one(A, j) ((A).coeff[(j) / BITS_PER_LONG] |= (1UL << ((j) % BITS_PER_LONG)))
#define D_mat_change_coeff(A, j) ((A).coeff[ (j) / BITS_PER_LONG] ^= (1UL << ((j) % BITS_PER_LONG)))
// Quasi_Dyadic Operation
#define QD_swaprows(A, i, k) {D_mat * pt = A.dblock[i]; A.dblock[i] =A.dblock[k]; A.dblock[k] = pt;}

binmat_t mat_ini(int rown, int coln);
D_mat D_mat_ini(int size);
QD_mat QD_mat_ini(int rown, int coln, int ordre);
QD_mat QD_mat_ini_Id(int rown, int coln, int ordre);
D_mat D_mat_ini_Id(int ordre);
void D_mat_aff(D_mat A, int ordre);
void QD_mat_aff(QD_mat A);
D_mat D_standard_multiplication_binary(int ordre, D_mat a, D_mat b);
D_mat D_karatsuba_binary(int ordre, D_mat a, D_mat b);
D_mat D_multiplication_binary(int ordre, D_mat a, D_mat b);
int D_binary_det(D_mat a, int ordre);
D_mat D_binary_sum(int ordre, D_mat a, D_mat b);
//int binary_quasidyadic_pivot(QD_mat a,int j, unsigned int* p);
int is_null_dyadic(D_mat a);
//QD_mat QD_LU_inverse(QD_mat M);
QD_mat multiply_quasi_dyadic_matrices(QD_mat A, QD_mat B);
//QD_mat QD_copy_sub_square(QD_mat M);
//int binary_quasidyadic_Square_LU(QD_mat M, unsigned int* permutation);
int D_GaussElim(QD_mat H);
binmat_t D_full_binary_from_Dyadic(QD_mat QD);
void QD_secret_matrix(QD_mat H, gf_t *u, gf_t *v, gf_t *z);
void QD_mat_pub(QD_mat R, QD_mat H);

//binmat_t mat_ini_from_string(int rown, int coln, const unsigned char * s);
void mat_free(binmat_t A);
void D_mat_free(D_mat A);
void QD_mat_free(QD_mat A);
//binmat_t mat_copy(binmat_t A);
void mat_rowxor(binmat_t A, int a, int b);
int* mat_rref(binmat_t A);
//void mat_vec_mul(unsigned long *cR, unsigned char *x, binmat_t A);
//binmat_t mat_mul(binmat_t A, binmat_t B);
void mat_aff(binmat_t A);
mat_t matrix_init(int rown, int coln);
void aff_mat(mat_t mat);

void secret_dyadic_bin_matrix(binmat_t H, gf_t *u, gf_t *v, gf_t *z);
void mat_free_nb(mat_t A);
void aff_bin_mat(binmat_t mat);
binmat_t mat_from_pk(const unsigned char* pk);

void produit_vector_matrix(unsigned long* Res, unsigned char* u, binmat_t A);
void G_mat_pub(binmat_t G, binmat_t H_syst);
void secret_matrix_new1(binmat_t H_bin, gf_t *u, gf_t *v, gf_t *z);
int gausselim(binmat_t M);
#endif

