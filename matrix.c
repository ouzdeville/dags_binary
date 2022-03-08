/*********************************************************************************************
 * BINARY DAGS: Key Encapsulation using binary Dyadic GS Codes.                            *
 * This is the binary version of DAGS  submitted to the NIST Post=Quantum Cryptography.    *
 * For any other usage , contact the author(s) to ascode_dimmension permission.            *
 **********************************************************************************************

 ~~~~~~~~Matrix Operations ~~~~~~~~~~~~~~~~ */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "matrix.h"
#include "time.h"

/*********************************************************************************************/
////////////////////////////////////MATRIX Functions///////////////////////////////////////////
/*********************************************************************************************/

//take a look at the MSB-LSB format???
binmat_t mat_ini(int rown, int coln) {
	binmat_t A;
	unsigned int i;
	//A = (binmat_t) malloc(sizeof(struct matrix));
	A.coln = coln;
	A.rown = rown;
	A.rwdcnt = (1 + (coln - 1) / BITS_PER_LONG);
	A.alloc_size = rown * A.rwdcnt * sizeof(unsigned long);
	A.elem = (unsigned long **) malloc(rown * sizeof(unsigned long *));
	for (i = 0; i < rown; i++)
		A.elem[i] = (unsigned long*) calloc(1 + (coln - 1) / BITS_PER_LONG,
				sizeof(unsigned long));
	return A;
}

mat_t matrix_init(int rown, int coln) {
	unsigned int i;
	mat_t A;
	A.coln = coln;
	A.rown = rown;

	A.coeff = (gf_t **) malloc(A.rown * sizeof(gf_t *));
	for (i = 0; i < A.rown; i++) {
		A.coeff[i] = (gf_t *) calloc(A.coln, sizeof(gf_t));
	}

	return A;
}

/**
 * TO initialize a signature for Dyadic Matrix
 */
D_mat D_mat_ini(int ordre) {
	D_mat a;
	//a.ordre=ordre;
	a.rwdcnt = (1 + (ordre - 1) / BITS_PER_LONG);
	a.alloc_size = a.rwdcnt * sizeof(unsigned long);
	a.coeff = (unsigned long*) calloc(1 + (ordre - 1) / BITS_PER_LONG,
			sizeof(unsigned long));
	return a;
}

/**
 * To initialize the signature of Ididentity as dyadic matrix
 */
D_mat D_mat_ini_Id(int ordre) {
	D_mat a;
	//a.ordre=ordre;
	a.rwdcnt = (1 + (ordre - 1) / BITS_PER_LONG);
	a.alloc_size = a.rwdcnt * sizeof(unsigned long);
	a.coeff = (unsigned long*) calloc(1 + (ordre - 1) / BITS_PER_LONG,
			sizeof(unsigned long));
	D_mat_set_coeff_to_one(a, 0);
	return a;
}

QD_mat QD_mat_ini(int rown, int coln, int ordre) {
	unsigned int i, j;
	QD_mat A;
	A.coln = coln;
	A.rown = rown;
	A.ordre = ordre;

	A.dblock = (D_mat **) malloc(rown * sizeof(D_mat *));
	for (i = 0; i < rown; i++) {
		A.dblock[i] = (D_mat *) calloc(coln, sizeof(D_mat));
		for (j = 0; j < coln; j++) {
			A.dblock[i][j] = D_mat_ini(order);
		}
	}

	return A;
}

QD_mat QD_mat_ini_Id(int rown, int coln, int ordre) {
	unsigned int i, j;
	QD_mat A;
	A.coln = coln;
	A.rown = rown;
	A.ordre = ordre;

	A.dblock = (D_mat **) malloc(rown * sizeof(D_mat *));
	for (i = 0; i < rown; i++) {
		A.dblock[i] = (D_mat *) calloc(coln, sizeof(D_mat));
		for (j = 0; j < coln; j++) {
			A.dblock[i][j] = D_mat_ini(order);
		}
		D_mat_set_coeff_to_one(A.dblock[i][i], 0);
	}

	return A;
}

void mat_free(binmat_t A) {
	unsigned int i;
	for (i = 0; i < A.rown; i++) {
		free(A.elem[i]);
	}
	free(A.elem);
}
void D_mat_free(D_mat A) {
	free(A.coeff);
}

void QD_mat_free(QD_mat A) {
	unsigned int i, j;
	for (i = 0; i < A.rown; i++) {
		for (j = 0; j < A.coln; ++j) {
			D_mat_free(A.dblock[i][j]);
		}

	}
	free(A.dblock);
}

void mat_rowxor(binmat_t A, int a, int b) {
	int i;
	for (i = 0; i < A.rwdcnt; i++) {
		A.elem[a][i] ^= A.elem[b][i];
	}
	//return A;
}

binmat_t mat_swapcol(binmat_t A, int a, int b) {
	int i;
	for (i = 0; i < A.rown; i++) {
		if (mat_coeff(A, i, a) != mat_coeff(A, i, b)) {
			mat_change_coeff(A, i, a);
			mat_change_coeff(A, i, b);
		}

	}
	return A;
}

//the matrix is reduced from LSB...(from right)

int* mat_rref(binmat_t A) {
	//printf("tour %d-%d\n",A->rown,A->coln);

	int i, j, findrow, first_col = A.coln - A.rown;
	int *perm;
	int col_ref = first_col;

	perm = malloc(A.coln * sizeof(int));

	for (i = 0; i < A.coln; i++)
		perm[i] = i; //initialize permutation.
	//failcnt = 0;

	for (i = 0; i < A.rown; i++) {
		findrow = 0;
		col_ref = first_col + i;

		for (j = i; j < A.rown; j++) {
			//printf("tour %d-%d",j,max);
			if (mat_coeff(A, j, col_ref))  //(A->elem[(j*A->coln)+max])
					{
				//max--;
				if (i != j)  //not needed as ith row is 0 and jth row is 1.
					swaprows(A, i, j);  //xor to the row.(swap)?
				findrow = 1;
				break;
			}  //largest value found (end if)
			   //		  break;
		}	  //fin j

		if (!findrow)//if no row with a 1 found then swap last column and the column with no 1 down.
		{

		}
		for (j = i + 1; j < A.rown; j++)	//fill the column downwards with 0's
				{
			if (mat_coeff(A, j, (col_ref)))	  //(A->elem[j*A->coln+max+1])
				mat_rowxor(A, j, i);	  //check the arg. order.
		}

		for (j = i - 1; j >= 0; j--)	 //fill the column with 0's upwards too.
				{
			if (mat_coeff(A, j, (col_ref)))	  //(A->elem[j*A->coln+max+1])
				mat_rowxor(A, j, i);
		}

	}	  //end for(i)

	return perm;
}

int gausselim(binmat_t M) {
	int i, j, k, l;
	//int rwdcnt = 1 + ( - 1) / __WORDSIZE;

	for (j = M.coln - M.rown, k = 0; (k < M.rown) && (j < M.coln); ++j) {
		// we search for a pivot
		for (i = k; i < M.rown; ++i) {
			if (mat_coeff(M, i, j)) {
				swaprows(M, i, k);

				break;
			}
		}

		if (i < M.rown) {      // found
			for (l = 0; l < k; ++l)
				if (mat_coeff(M, l, j))
					mat_rowxor(M, l, k);
			for (l = i + 1; l < M.rown; ++l)
				if (mat_coeff(M, l, j))
					mat_rowxor(M, l, k);
			++k;
		} else {          // not found

		}
	}

	return k;
}

int gaussmc(binmat_t M) {

	int i, j, k;
	int row, c;
	unsigned long mask;

	for (i = 0; i < ((M.rown - 1) / BITS_PER_LONG) + 1; i++)
		for (j = 0; j < BITS_PER_LONG; j++) {
			row = i * BITS_PER_LONG + j;

			if (row >= M.rown)
				break;

			for (k = row + 1; k < M.rown; k++) {
				mask = M.elem[row][i] ^ M.elem[k][i];
				mask >>= j;
				mask &= 1;
				mask = -mask;

				for (c = 0; c < M.rwdcnt; c++)
					M.elem[row][c] ^= M.elem[k][c] & mask;
			}

			if (((M.elem[row][i] >> j) & 1) == 0) // return if not systematic
					{
				return -1;
			}

			for (k = 0; k < M.rown; k++) {
				if (k != row) {
					mask = M.elem[k][i] >> j;
					mask &= 1;
					mask = -mask;

					for (c = 0; c < M.rwdcnt; c++)
						M.elem[k][c] ^= M.elem[row][c] & mask;
				}
			}
		}
	return 0;
}

void mat_aff(binmat_t A) {
	int i, k;
	for (i = 0; i < A.rown; i++) {
		for (k = 0; k < A.coln; k++)

		{
			//printf("%d\t", mat.coeff[i][j]);
			printf("%lu", mat_coeff(A, i, k));

		}
		printf("\n");
	}

	printf("\n");
}

void aff_mat(mat_t mat) {
	printf("Show new Matrix\n");
	int i, j;
	for (i = 0; i < (mat.rown); i++) {
		for (j = 0; j < (mat.coln); j++) {
			printf("%d-", mat.coeff[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

void secret_dyadic_bin_matrix(binmat_t H, gf_t *u, gf_t *v, gf_t *z) {

	gf_t *Z;
	Z = (gf_t *) calloc(code_length, sizeof(gf_t));
	int i, j, k, l = 0;
	gf_t hij = 1, h = 1;
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	///// Zi fillfull the following restriction: Zis+j=Zis, for i=0...n0-1, j=0...s-1
	//// That means that Z0=Z1=...=Z(s-1); Zs=Z(s+1)=...=Z(2s-1);  ... Z(n-s-1)=Z(n-s)=...=Z(n-1)
	//// At the end, we just have to choose n0 distincts elements.
	///////////////////////////////////////////////////////////////////////////////////////////////////////

	for (i = 0; i < n0_val; i++) {
		for (j = 0; j < order; j++) {
			Z[i * (order) + j] = z[i];
		}
	}

	for (j = 0; j < code_length; ++j) {

		for (i = 0; i < order; ++i) {
			hij = gf_inv(u[i] ^ v[j]);
			h = Z[j];
			for (l = 0; l < pol_deg; ++l) {
				h = gf_mul(h, hij);
				for (k = 0; k < EXT_DEGREE; ++k) {
					if (h & (1 << k)) {
						mat_set_coeff_to_one(H, k*pol_deg*order+l*order+i, j); //the co-eff. are set in 2^0,...,2^11 ; 2^0,...,2^11 format along the rows/cols?
					}
				}
			}

		}
	}
}

void mat_free_nb(mat_t A) {
	unsigned int i;
	for (i = 0; i < A.rown; i++) {
		free(A.coeff[i]);
	}
	free(A.coeff);
}
void aff_bin_mat(binmat_t mat) {
	printf("Show new Bin Matrix\n");
	int i, j;
	for (i = 0; i < (mat.rown); i++) {
		for (j = 0; j < (mat.coln); j++) {
			printf("%lu", mat_coeff(mat, i, j));
		}
		printf("\n");
	}
	printf("\n");
}

binmat_t mat_from_pk(const unsigned char* pk) {
	binmat_t H;
	int r = EXT_DEGREE * pol_deg * order;
	int n = code_length;
	H = mat_ini(n - r, r);

	int i, j, k;
	int ind;
	for (i = 0; i < H.rown; i += order) {
		memcpy(H.elem[i], pk, H.rwdcnt * sizeof(unsigned long));
		pk += H.rwdcnt * sizeof(unsigned long);
		for (k = i + 1; k < i + order; k++) {
			for (j = 0; j < r; j++) {
				ind = k ^ j ^ i;
				if (mat_coeff(H, i, ind)) {
					mat_set_coeff_to_one(H, k, j);

				}

			}

		}

	}
	return H;
}

void produit_vector_matrix(unsigned long* Res, unsigned char* u, binmat_t A) {
	int i, j;
	for (i = 0; i < A.rown; i++) {
		if ((u[i / 8] >> (i % 8)) & 1) {
			for (j = 0; j < A.rwdcnt; j++) {
				Res[j] ^= A.elem[i][j];
			}
		}

	}

}

void G_mat_pub(binmat_t R, binmat_t H_syst) {
	int i, k;
	for (i = 0; i < R.rown; i++) {
		for (k = 0; k < R.coln; k++)

		{
			//printf("%d\t", mat.coeff[i][j]);
			if (mat_coeff(H_syst, k, i))
				mat_set_coeff_to_one(R, i, k);

		}
	}

}

void secret_matrix_new1(binmat_t H_bin, gf_t *u, gf_t *v, gf_t *z) {
	mat_t T[pol_deg]; // T will contain all the block matrices H1 to H2
	gf_t *Z;
	Z = (gf_t *) calloc(code_length, sizeof(gf_t));
	int i, j, k, s = order, t = pol_deg, l = 0;
	gf_t temp = 0;
	//mat_set_to_zero(H_bin);

	for (i = 0; i < n0_val; i++) {
		for (j = 0; j < order; j++) {
			Z[i * (order) + j] = z[i];
		}
	}

	for (i = 0; i < pol_deg; i++) {
		T[i] = matrix_init(order, code_length);
	}

	for (i = 0; i < order; i++) {
		for (j = 0; j < code_length; j++) {
			T[0].coeff[i][j] = gf_inv(u[i] ^ v[j]);
			temp = gf_mul(T[0].coeff[i][j], Z[j]);
			for (k = 0; k < gf_extd(); k++) {
				if (temp & (1 << k)) {
					mat_set_coeff_to_one(H_bin, k * t * (s) + i, j);

				}
			}
		}
	}
	for (j = 0; j < code_length; j++) {
		for (i = 0; i < order; i++) {
			for (l = 1; l < pol_deg; l++) {
				T[l].coeff[i][j] = gf_mul(gf_inv(u[i] ^ v[j]),
						T[l - 1].coeff[i][j]);
				temp = gf_mul(T[l].coeff[i][j], Z[j]);
				for (k = 0; k < gf_extd(); k++) {
					if (temp & (1 << k)) {
						mat_set_coeff_to_one(H_bin, k * t * (s) + l * s + i, j);

					}

				}

			}
		}
	}

	free(Z);
	for (i = 0; i < pol_deg; i++) {
		mat_free_nb(T[i]);
	}
}

void D_mat_aff(D_mat A, int ordre) {
	int k;
	for (k = 0; k < ordre; k++)

	{
		//printf("%d\t", mat.coeff[i][j]);
		printf("%lu", D_mat_coeff(A, k));

	}
	printf(" ");

}
void QD_mat_aff(QD_mat mat) {
	printf("Show new QD_Matrix\n");
	int i, j;
	for (i = 0; i < (mat.rown); i++) {
		for (j = 0; j < (mat.coln); j++) {
			D_mat_aff(mat.dblock[i][j], mat.ordre);
		}
		printf("\n\n");
	}
	printf("\n");
}

void QD_secret_matrix(QD_mat H, gf_t *u, gf_t *v, gf_t *z) {
	gf_t *Z;
	Z = (gf_t *) calloc(code_length, sizeof(gf_t));
	int i, j, k, l = 0;
	gf_t hij = 1, h = 1;
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	///// Zi fillfull the following restriction: Zis+j=Zis, for i=0...n0-1, j=0...s-1
	//// That means that Z0=Z1=...=Z(s-1); Zs=Z(s+1)=...=Z(2s-1);  ... Z(n-s-1)=Z(n-s)=...=Z(n-1)
	//// At the end, we just have to choose n0 distincts elements.
	///////////////////////////////////////////////////////////////////////////////////////////////////////

	for (i = 0; i < n0_val; i++) {
		for (j = 0; j < order; j++) {
			Z[i * (order) + j] = z[i];
		}
	}

	for (j = 0; j < code_length; ++j) {

		for (i = 0; i < order; i = i + order) {
			hij = gf_inv(u[i] ^ v[j]);
			h = Z[j];
			for (l = 0; l < pol_deg; ++l) {
				h = gf_mul(h, hij);
				for (k = 0; k < EXT_DEGREE; ++k) {
					if (h & (1 << k)) {
						//mat_set_coeff_to_one(H, k*pol_deg*order+l*order+i, j); //the co-eff. are set in 2^0,...,2^11 ; 2^0,...,2^11 format along the rows/cols?
						//k*pol_deg*order+l*order+i
						D_mat_set_coeff_to_one(H.dblock[k*pol_deg+l][j/order],
								j%order);
						//D_mat_set_coeff_to_one(QD.dblock[i/order][j/order], j%order);
					}
				}
			}

		}
	}
}

D_mat D_standard_multiplication_binary(int ordre, D_mat a, D_mat b) {
	unsigned int i, j;
	D_mat c_sig = D_mat_ini(ordre);
	//D_mat_aff(c_sig,ordre);
	unsigned int A, B;
	for (i = 0; i < ordre; i++) {

		for (j = 0; j < ordre; ++j) {
			B = D_mat_coeff(b, i ^ j);
			A = D_mat_coeff(a, j);
			if ((A & B))
				D_mat_change_coeff(c_sig, i);
		}

	}
	return c_sig;
}

D_mat D_multiplication_binary(int ordre, D_mat a, D_mat b) {
	return D_standard_multiplication_binary(ordre, a, b);
}

/**
 * Karatsuba Multiplication of two binary matrices by its signatures
 * return Delta(a)*Delta(b)
 *
 */
D_mat D_karatsuba_binary(int ordre, D_mat a, D_mat b) {
	int i;
	int k = ordre;
	D_mat c_sig = D_mat_ini(k);
	//printf("Ordre=%d\n",ordre);
	if (k <= 2) {
		//c_sig.[0] = gf_mult_fast(a[0] , b[0]) ^ gf_mult_fast(a[1] , b[1]);
		//c_sig[1] = gf_mult_fast(a[0] , b[1]) ^ gf_mult_fast(a[1] , b[0]);
		if ((D_mat_coeff(a, 0) & D_mat_coeff(b, 0))
				^ (D_mat_coeff(a, 1) & D_mat_coeff(b, 1)))
			D_mat_set_coeff_to_one(c_sig, 0);
		if ((D_mat_coeff(a, 0) & D_mat_coeff(b, 1))
				^ (D_mat_coeff(a, 1) & D_mat_coeff(b, 0)))
			D_mat_set_coeff_to_one(c_sig, 1);
		return c_sig;

	} else {
		int k_2 = k / 2;
		//First half of the result
		D_mat a0 = D_mat_ini(k_2);
		D_mat b0 = D_mat_ini(k_2);
		D_mat c0 = D_mat_ini(k_2);
		D_mat a1 = D_mat_ini(k_2);
		D_mat b1 = D_mat_ini(k_2);
		D_mat c1 = D_mat_ini(k_2);

		for (int i = 0; i < k_2; i++) {
			if (D_mat_coeff(c_sig, i))
				D_mat_set_coeff_to_one(c0, i);

			if (D_mat_coeff(a, i))
				D_mat_set_coeff_to_one(a0, i);

			if (D_mat_coeff(b, i))
				D_mat_set_coeff_to_one(b0, i);

			if (D_mat_coeff(c_sig, i + k_2))
				D_mat_set_coeff_to_one(c1, i);

			if (D_mat_coeff(a, i + k_2))
				D_mat_set_coeff_to_one(a1, i);

			if (D_mat_coeff(b, i + k_2))
				D_mat_set_coeff_to_one(b1, i);
		}
		//c0=a0*b0
		c0 = D_karatsuba_binary(k_2, a0, b0);

		//c1=a1*b1
		c1 = D_karatsuba_binary(k_2, a1, b1);

		//a0.b0+a1.b1 // to be optimized
		for (i = k_2; i--;) {
			if ((D_mat_coeff(c0, i)) ^ (D_mat_coeff(c1, i)))
				D_mat_set_coeff_to_one(c_sig, i);
		}

		D_mat a0_plus_a1 = D_mat_ini(k_2);
		D_mat b0_plus_b1 = D_mat_ini(k_2);
		D_mat a01_prod_b01 = D_mat_ini(k_2);

		a0_plus_a1 = D_binary_sum(k_2, a0, a1);
		b0_plus_b1 = D_binary_sum(k_2, b0, b1);

		a01_prod_b01 = D_karatsuba_binary(k_2, a0_plus_a1, b0_plus_b1);

		for (int i = k_2; i < k; i++) {
			if ((D_mat_coeff(a01_prod_b01, i - k_2))
					^ (D_mat_coeff(c_sig, i - k_2)))
				D_mat_set_coeff_to_one(c_sig, i);
		}
	}
	return c_sig;

}

/*
 * Computation of de determinant of binary dyadic matrix  by its signature
 */
int D_binary_det(D_mat a, int ordre) {
	int i;
	unsigned int det = 0;

	for (i = 0; i < ordre; i++) {
		det = det ^ D_mat_coeff(a, i);
	}
	return det;
}

/**
 * Sum of two dyadic binary matrices by its signatures
 */
D_mat D_binary_sum(int ordre, D_mat a, D_mat b) {
	int i;
	D_mat c = D_mat_ini(ordre);
	for (i = 0; i < a.rwdcnt; i++) {
		c.coeff[i] = a.coeff[i] ^ b.coeff[i];
	}
	return c;
}

/**
 * tells if a binary dyadic matrix is not null?
 */
int is_not_null_dyadic(D_mat a) {
	int i;
	int val = 0;
	for (i = 0; i < a.rwdcnt; i++) {
		val += a.coeff[i];
	}
	return val;
}
QD_mat multiply_quasi_dyadic_matrices(QD_mat A, QD_mat B) {
	QD_mat res;
	if (A.coln != B.rown) {
		printf("Error: Impossible to multiply\n");
		return res;
	} else {
		int i, j, k;

		D_mat a;
		res = QD_mat_ini(A.rown, B.coln, A.ordre);
		for (i = 0; i < A.rown; i++) {
			for (j = 0; j < B.coln; j++) {
				a = D_mat_ini(A.ordre);
				for (k = 0; k < A.coln; k++) {
					a = D_binary_sum(A.ordre, a,
							D_multiplication_binary(A.ordre, A.dblock[i][k],
									B.dblock[k][j]));
				}
				res.dblock[i][j] = a;
			}
		}
		return res;
	}
}

int D_GaussElim(QD_mat H) {

	int i, j, l = 0, test = 0;
	D_mat* temp;
	int n = H.coln;
	int k = H.rown;
	int w, nb = n - k;
	for (i = 0; i < k; i++) {

		//QD_mat_aff(H);
		test = 0;

		l = 0;
		w = i + nb;
		j = w;
		//printf("I,J= %d,%d\n", i,j);
		if (!D_binary_det(H.dblock[i][w], H.ordre)) { //We're looking for a non-invertible pivot
			test = 1;
			//printf("search Pivot\n");
			for (l = i + 1; l < k; l++) {
				//printf("L= %d\n", l);
				if (D_binary_det(H.dblock[l][j], H.ordre)) {
					//printf("Find Pivot\n");
					break;
				}
			}
		}

		if (l == k) {

			printf("Non systematic D-Matrix %d\n", l);

			return -1;
		}
		if (test == 1) { // We switches the lines l and i
			test = 0;
			//printf("Permut line\n");
			//temp=P[i+n-k];
			//P[i+n-k]=P[j];
			//P[j]=temp;
			/*for (j = 0; j < n; j++)
			 {
			 temp = H.coeff[l][j];
			 H.coeff[l][j] = H.coeff[i][j];
			 H.coeff[i][j] = temp;
			 }*/
			//printf("AfTER Pivot foud l=%d\n",l);
			temp = H.dblock[l];
			H.dblock[l] = H.dblock[i];
			H.dblock[i] = temp;
		}
		//   Matrix standardization
		D_mat invPiv = D_mat_ini_Id(H.ordre);
		//D_mat_aff(invPiv,H.ordre);
		//printf("INV\n");
		//si la matrice n'est pas l'identité Pivot on normalise
		if (is_not_null_dyadic(D_binary_sum(H.ordre, H.dblock[i][w], invPiv))) {
			//D_mat_aff(H.dblock[i][w],H.ordre);
			//printf("\n");
			//D_mat_aff(invPiv,H.ordre);

			for (j = 0; j < n; j++) {
				//printf("inv=");
				//D_mat_aff(H.dblock[i][j],H.ordre);
				//D_mat_aff(H.dblock[i][w],H.ordre);
				H.dblock[i][j] = D_multiplication_binary(H.ordre,
						H.dblock[i][j], H.dblock[i][w]);
				//D_mat_aff(H.dblock[i][j],H.ordre);
				//printf("\n");
			}

		}

		//Here we do the elimination on column i + n-k
		D_mat piv_align;
		for (l = 0; l < k; l++) {
			if (l == i) {
				continue;
			}
			piv_align = H.dblock[l][w];
			if (is_not_null_dyadic(piv_align)) {

				for (j = 0; j < n; j++) {
					H.dblock[l][j] = D_binary_sum(H.ordre, H.dblock[l][j],
							(D_multiplication_binary(H.ordre, piv_align,
									H.dblock[i][j])));
					// piv_align is useless because it does before, so you can use gf_mul_fast_subfield directly
				}
			}
		}
	}

	return H.rown;
}

binmat_t D_full_binary_from_Dyadic(QD_mat QD) {
	int i, j;
	binmat_t Hs = mat_ini(QD.rown, QD.coln * order); //initialize matrix with actual no. of bits.
	//printf("In Secret PCM in QD FORM\n");
	for (i = 0; i < Hs.rown; i++) {
		for (j = 0; j < Hs.coln; ++j) {
			if (D_mat_coeff(QD.dblock[i][j/order], j%order)) {
				mat_set_coeff_to_one(Hs, i, j);
			}
		}
	}
	return Hs;
}
void QD_mat_pub(QD_mat R, QD_mat H) {
	int i, k, j;
	for (i = 0; i < R.rown; i++) {
		for (k = 0; k < R.coln; k++)

		{
			//printf("%d\t", mat.coeff[i][j]);
			if (is_not_null_dyadic(H.dblock[k][i])) {
				for (j = 0; j < H.dblock[k][i].rwdcnt; ++j) {
					R.dblock[i][k].coeff[j] = H.dblock[k][i].coeff[j];
				}
			}

		}
	}

}

