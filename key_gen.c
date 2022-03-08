/******************************************************************************************
 * BINARY DAGS: Key Encapsulation using binary Dyadic GS Codes.                            *
 * This is the binary version of DAGS  submitted to the NIST Post=Quantum Cryptography.    *
 * For any other usage , contact the author(s) to ascode_dimmension permission.            *
 *******************************************************************************************
 */

#include "key_gen.h"

int disjoint_test(gf_t *u, gf_t *v) {
	int i, j;
	for (i = 0; i < (order); i++) {
		for (j = 0; j < code_length; j++) {
			if (u[i] == v[j]) {
				return -1;
			}
		}
	}
	return 0;
}

int Test_disjoint(gf_t *L, int n) {
	int i, j;
	for (i = 0; i < n; i++) {
		for (j = i + 1; j < n; j++) {
			if (L[i] == L[j]) {
				return -1;
			}
		}
	}
	return 0;
}

void generate_random_vector(int m, gf_t *vect) {
	int i, j, v;
	gf_t tmp;
	gf_t *U;
	U = (gf_t *) calloc(gf_card(), sizeof(gf_t));
	unsigned char *random_bytes = malloc(gf_card() * sizeof(gf_t));
	randombytes(random_bytes, gf_card() * sizeof(gf_t));
	U[0] = 1;
	for (i = 1; i < gf_card(); i++) {
		U[i] = i;
	}
	for (j = 1; j < gf_card(); j++) {

		v = ((gf_t *) random_bytes)[j] % (j + 1);
		tmp = U[j];
		U[j] = U[v + 1];
		U[v + 1] = tmp;
	}

	memcpy(vect, U + 1, (m) * sizeof(gf_t));

	free(U);
	free(random_bytes);
}

void init_random_element(gf_t *U) {
	int i, j, v;
	gf_t tmp;
	unsigned char *random_bytes = 0;
	random_bytes = malloc(gf_ord() * sizeof(gf_t));
	randombytes(random_bytes, gf_ord() * sizeof(gf_t));
	for (i = 0; i <= gf_ord(); i++) {
		U[i] = i;
	}

	for (j = 1; j < gf_ord(); j++) {

		v = ((gf_t *) random_bytes)[j] % (j + 1);
		tmp = U[j];
		U[j] = U[v + 1];
		U[v + 1] = tmp;
	}

	free(random_bytes);
}

void Remove_From_U(gf_t elt, gf_t *U) {
	int k;
	for (k = 0; k <= gf_ord(); k++) {
		if (U[k] == elt) {
			U[k] = 0;
			break;
		}
	}
}

void binary_quasi_dyadic_sig(int m, int n, int t, int *b, gf_t *h_sig, gf_t *w) {
	int i, j, k, s, p, l, l1, c, r, consistent_root, consistent_support_block;
	int const C = ((gf_card()) / t);
	gf_t *U, *V, *h;
	gf_t sum_inv_h_i_j_0, sum_inv_h_i_0;
	h = (gf_t *) calloc(gf_card(), sizeof(gf_t));
	U = (gf_t *) calloc(gf_card(), sizeof(gf_t));
	V = (gf_t *) calloc(gf_card(), sizeof(gf_t));

	do {
		init_random_element(U);
		h[0] = U[1];
		U[1] = 0;

		for (s = 0; s < m; s++) {
			i = 1 << s;
			h[i] = U[i + 1];
			Remove_From_U(h[i], U);
			for (j = 1; j < i; j++) {
				h[i + j] = 0;
				if ((h[i] != 0) && (h[j] != 0)) {
					sum_inv_h_i_j_0 = (gf_inv(h[i])^ gf_inv(h[j])) ^ (gf_inv(h[0]));
					if (sum_inv_h_i_j_0 != 0)
					{
						h[i + j] = gf_inv(sum_inv_h_i_j_0);
						Remove_From_U(h[i + j], U);
					}
					else
					{
						h[i + j] = 0;
					}
				}
				else
				{
					h[i + j] = 0;
				}
			}
		}

		c = 0;
		init_random_element(V);
		consistent_root = 1;
		for (p = 0; p < t; p++) {
			consistent_root = consistent_root & (h[p] != 0);
		}
		if (consistent_root) {
			b[0] = 0;
			c = 1;
			for (r = 0; r < t; r++) {
				sum_inv_h_i_0 = (gf_inv(h[r])) ^ (gf_inv(h[0]));
				Remove_From_U(gf_inv(h[r]), V);
				Remove_From_U(sum_inv_h_i_0, V);
			}
			for (j = 1; j < C; j++) {
				consistent_support_block = 1;
				for (p = j * t; p < (j + 1) * t; p++) {
					consistent_support_block = consistent_support_block
							& (h[p] != 0);
				}
				if (consistent_support_block) {
					b[c] = j;
					c = c + 1;
					for (l = j * t; l < (j + 1) * t; l++) {
						sum_inv_h_i_0 = (gf_inv(h[l])) ^ (gf_inv(h[0]));
						Remove_From_U(sum_inv_h_i_0, V);
					}
				}
			}
		}
	} while (c * t < n);

	// Computing w: We just one value of omega. So we stop at the first non-zero element of V.
	for (j = 0; j < gf_card(); j++) {
		if (V[j]) {
			*w = V[j];
			break;
		}
	}
	/******************************************
	 We choose n0=33 consistent blocks from all the consistent blocks given by the vector  b;
	 We then obtain
	 ******************************************/
	for (j = 0; j < n0_val; j++) {
		for (k = 0; k < order; k++) {
			l = (order) * j + k;
			l1 = (order) * b[j] + k;
			h_sig[l] = h[l1];
		}
	}

	free(U);
	free(V);
	free(h);
}

void cauchy_support(gf_t *Support, gf_t *W, gf_t *w) {
	int i;
	gf_t sum_inv_h_i_0;
	gf_t *h;
	int *b, test_u = 0, test_v = 0, test_u_inter_v = 0;
	do {
		b = (int *) calloc(gf_card(), sizeof(int));
		h = (gf_t *) calloc(code_length, sizeof(gf_t));
		binary_quasi_dyadic_sig(gf_extd(), code_length, order, b, h, w);
		for (i = 0; i < code_length; i++) {
			sum_inv_h_i_0 = (gf_inv(h[i])) ^ (gf_inv(h[0]));
			Support[i] = (sum_inv_h_i_0) ^ (w[0]);
		}
		for (i = 0; i < order; i++) {
			W[i] = (gf_inv(h[i])) ^ (w[0]);
		}
		test_u = Test_disjoint(Support, code_length);
		test_v = Test_disjoint(W, order);
		test_u_inter_v = disjoint_test(W, Support);

	} while ((test_u != 0) || (test_v != 0) || (test_u_inter_v != 0));

	free(h);
	free(b);
}
/**
 *
 */

/*
 * The function key_pair generates the public key and the
 * secret key which will be stored in files
 */
int key_pair(unsigned char *pk, unsigned char *sk) {

	gf_t *u, *v, *w, *z;
	int return_value = 1;
	binmat_t Hs;
	binmat_t R;

	gf_init(EXT_DEGREE);

	int i;
	int n = code_length;
	int r = pol_deg * (order) * EXT_DEGREE;
	Hs = mat_ini(r, n); //initialize matrix with actual no. of bits.
	//QD_mat QD;

	gf_t* y = (gf_t *) calloc(code_length, sizeof(gf_t));
	while (return_value != Hs.rown) {
		u = (gf_t *) calloc(order, sizeof(gf_t));
		v = (gf_t *) calloc(code_length, sizeof(gf_t));
		w = (gf_t *) calloc(code_length, sizeof(gf_t));
		z = (gf_t *) calloc(n0_val, sizeof(gf_t));

		cauchy_support(v, u, w);
		//mat_aff(Hs);
		free(w);
		generate_random_vector(n0_val, z);

		set_y_from_uvz(u, v, z, y);
		set_sk_from_vy(v, y, sk);
		Hs = mat_ini(r, n);
		//QD = QD_mat_ini(mst/order, n0_val, order);
		secret_dyadic_bin_matrix(Hs, u, v, z);
		//QD_secret_matrix(QD, u, v, z);
		//QD_mat_aff(QD);
		//mat_aff(Hs);

		return_value = gausselim(Hs);
		//return_value = D_GaussElim(QD);

		if (return_value != Hs.rown) {
			mat_free(Hs);
			//QD_mat_free(QD);
		}

	}

	//mat_aff(Hs);
	//QD_mat_aff(QD);

	//Dealing Quasi-dyadic Block
	/*QD_mat QD_R = QD_mat_ini(QD.coln - QD.rown, QD.rown, QD.ordre);
	 QD_mat_pub(QD_R,QD);
	 QD_mat_free(QD);
	 //QD_mat_aff(QD_R);

	 binmat_t R = D_full_binary_from_Dyadic(QD_R);
	 //mat_aff(R);
	 for (i = 0; i < R.rown; i ++) {
	 memcpy(pk, R.elem[i], R.rwdcnt * sizeof(unsigned long));
	 pk += R.rwdcnt * sizeof(unsigned long);
	 //printf("%lu+", i);
	 }
	 pk -= R.rown * R.rwdcnt * sizeof(unsigned long);

	 mat_free(R);*/

	//Dealing with full binary matrix
	R = mat_ini(n - r, r);

	G_mat_pub(R, Hs);
	mat_free(Hs);

	for (i = 0; i < R.rown; i += order) {
		memcpy(pk, R.elem[i], R.rwdcnt * sizeof(unsigned long));
		pk += R.rwdcnt * sizeof(unsigned long);

	}

	mat_free(R);

	return 0;
}
