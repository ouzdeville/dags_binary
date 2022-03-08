/******************************************************************************************
 * BINARY DAGS: Key Encapsulation using binary Dyadic GS Codes.                            *
 * This is the binary version of DAGS  submitted to the NIST Post=Quantum Cryptography.    *
 * For any other usage , contact the author(s) to ascode_dimmension permission.            *
 *******************************************************************************************
 */

#include "encapsulation.h"
#include "api.h"

int encapsulation(const unsigned char *pk, unsigned char *ct, unsigned char *ss) {
	gf_init(EXT_DEGREE);
	int i;

	/*Step_1:  Choose randomly  m ←  F_q^k, m is seen as a sequence of k_prime integer modulo 2^6.......
	 *****************************************************************************************************/
	//variable declaration
	unsigned char* m = (unsigned char*) calloc(BITS_TO_BYTES(k_prime),
			sizeof(unsigned char));

	const unsigned char *custom = (unsigned char *) "DAGs"; // customization = "DAGs";

	//choose m randomly by using  by using random_m function.

	m = random_m(BITS_TO_BYTES(k_prime));

	/*Step_2:  Compute r = G(m) and d = H(m) with  G(x) = sponge(x,k) and H(x) = sponge(x,k_prime).........
	 *******************************************************************************************************/
	//variable declaration
	//gf* dd=( gf*)calloc(k_prime,sizeof( gf));
	unsigned char* r = (unsigned char*) calloc(BITS_TO_BYTES(code_dimension),
			sizeof(unsigned char));
	unsigned char* d = (unsigned char*) calloc(BITS_TO_BYTES(k_prime),
			sizeof(unsigned char));
	//r = (unsigned char *)calloc(BITS_TO_BYTES(code_dimension), sizeof(unsigned char));

	int test = KangarooTwelve(m, BITS_TO_BYTES(k_prime), r,
			BITS_TO_BYTES(code_dimension), custom, cus_len);
	assert(test == 0); // Catch Error

	test = KangarooTwelve(m, BITS_TO_BYTES(k_prime), d, BITS_TO_BYTES(k_prime),
			custom, cus_len);
	assert(test == 0); // Catch Error

	//for (i = 0; i < k_prime; i++)
	// Optimize modulo

	//dd[i] = (unsigned char)(d[i] & 1);

	//free(d);

	printf("\n ");

	/*Step_3:  Parse r as (ρ||σ) then set u = (ρ||m).......................................................
	 *******************************************************************************************************/
	//variable declaration
	unsigned char* rho = (unsigned char*) calloc(BITS_TO_BYTES(k_sec),
			sizeof(unsigned char));
	unsigned char* sigma = (unsigned char*) calloc(BITS_TO_BYTES(k_prime),
			sizeof(unsigned char));
	unsigned char* u = (unsigned char*) calloc(BITS_TO_BYTES(code_dimension),
			sizeof(unsigned char));

	/*
	 * Step_3:  Parse r as (ρ||σ) then set u = (ρ||m)
	 */
	memcpy(rho, r, BITS_TO_BYTES(k_sec));
	memcpy(sigma, r + BITS_TO_BYTES(k_sec), BITS_TO_BYTES(k_prime));

	memcpy(u, rho, BITS_TO_BYTES(k_sec));
	memcpy(u + BITS_TO_BYTES(k_sec), m, BITS_TO_BYTES(k_prime));

	free(r);
	free(rho);

	/*Step_4: Generate error vector e of length n and weight w from sigma.........................................
	 *********************************************************************************************************/
	//variable declaration
	unsigned char* hash_sigma = (unsigned char*) calloc(
			BITS_TO_BYTES(code_length), sizeof(unsigned char));

	// sigma: input type unsigned char len k_prime | hash_sigma: output type unsigned char len code_length
	test = KangarooTwelve(sigma, BITS_TO_BYTES(k_prime), hash_sigma,
			BITS_TO_BYTES(code_length), custom, cus_len);
	assert(test == 0); // Catch Error

	//Generate error vector e of length code_length and weight n0_w from hash_sigma1 by using random_e function.

	unsigned char *error_array;
	error_array = random_e(code_length, 2, n0_w, hash_sigma);

	free(sigma);
	free(hash_sigma);

	/*Step_5: Recovery of G and Compute c = uG + e................................................................
	 *************************************************************************************************************/
	//variable declaration
	unsigned long *c1 = (unsigned long*) calloc(
			1 + (code_length - code_dimension - 1) / BITS_PER_LONG,
			sizeof(unsigned long));
	//unsigned char* u1=( unsigned char*)calloc(BITS_TO_BYTES(code_dimension),sizeof( unsigned char));
	unsigned char* c = (unsigned char*) calloc(BITS_TO_BYTES(code_length),
			sizeof(unsigned char));

	// Construction of public matric from pk by using recup_pk function.

	int n = code_length;
	int r1 = pol_deg * (order) * EXT_DEGREE;
	binmat_t R = mat_ini(n - r1, r1);
	R = mat_from_pk(pk);

	// Compute code-word c with public matrix G and message u.
	//for(i = 0;i<code_dimension;i++) u1[i]=(gf)u[i];

	produit_vector_matrix(c1, u, R);
	memcpy(c, u, BITS_TO_BYTES(code_dimension));
	memcpy(c + BITS_TO_BYTES(code_dimension), c1,
			R.rwdcnt * sizeof(unsigned long));

	for (i = 0; i < code_length / 8; i++) {

		ct[i] = (unsigned char) ((unsigned char) c[i]
				^ (unsigned char) error_array[i]);

	}
	memcpy(ct + BITS_TO_BYTES(code_length), d, BITS_TO_BYTES(k_prime));

	free(u);
	free(c);
	// free(dd);
	free(error_array);

	/*Step_6: Compute K = K(m)...................................................................................
	 ************************************************************************************************************/

	//unsigned char *K = (unsigned char *)calloc(ss_length, sizeof(unsigned char));
	// Replace by KangarooTwelve
	// m: input type unsigned char len k_prime | K: output type unsigned char len ss_length
	test = KangarooTwelve(m, BITS_TO_BYTES(k_prime), ss, ss_length, custom,
			cus_len);
	assert(test == 0); // Catch Error

	printf(
			"\n*************************Encapsulation::La clef de session est:****************************\n\n");
	for (i = 0; i < ss_length; i++)
		printf("%d ", ss[i]);

	printf("\nFIN\n");

	//free(K);
	free(m);

	return 0;
	/*END********************************************************************************************************/

}

