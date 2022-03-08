/******************************************************************************************
* BINARY DAGS: Key Encapsulation using binary Dyadic GS Codes.                            *
* This is the binary version of DAGS  submitted to the NIST Post=Quantum Cryptography.    *
* For any other usage , contact the author(s) to ascode_dimmension permission.            *
*******************************************************************************************
*/

#include "decapsulation.h"

/*
 * Decapsulation() fuction compute the shared secret (ss) of type
 * unsigned char* and the ciphertext *(ct) of type unsigned char* by using the
 * secret key (sk)                             *
 */
int decapsulation(unsigned char *ss, const unsigned char *ct,
                  const unsigned char *sk)
{

    int i, test, decode_value;
    gf_init(EXT_DEGREE);                                            // Initialize of Log Antilog table
    const unsigned char *custom = (unsigned char *)"DAGs"; // customization = "DAGs";
    unsigned char *mot;
    unsigned char *m1, *rho1;
    unsigned char *r1, *d1, *rho2, *sigma, *e2, *hash_sigma, *e_prime;
    

        gf* v = (gf *)calloc(code_length, sizeof(gf));
        gf* y = (gf *)calloc(code_length, sizeof(gf));
	//gf* cipher = (gf *)calloc(code_length, sizeof(gf));
	
	//memcpy(v, sk, code_length*sizeof(gf_t));
        //memcpy(y,sk+code_length*sizeof(gf_t), sizeof(gf_t)*code_length);

	//unsigned char *ciph = ( unsigned char *)calloc(code_length, sizeof(unsigned char));

	


	printf("\n");

        set_vy_from_sk(v,y,sk);
	
	
    /*
    * Read in the alternative matrix from the secret key
    */
    

    /*
   * Step_1 of the decapsulation :  Decode the noisy codeword C received as
   * part of the ciphertext ct = (c||d) with d is “ a plaintext confirmation”.
   * We obtain codeword mot = u1G and error e
   */
    e_prime = (unsigned char *)calloc(BITS_TO_BYTES(code_length), sizeof(unsigned char));
    mot = (unsigned char *)malloc(BITS_TO_BYTES(code_length));

	printf("\nAFTER_DECODE\n");
    memcpy(mot, ct, BITS_TO_BYTES(code_length));
    decode_value =decoding_from_vy(v,y, ct, e_prime, mot);
    free(y);
    free(v);


    /*
   * Step_2 of the decapsulation :  Output ⊥ if decoding fails or wt(e) != n0_w
   */


    if (decode_value == -1 || bin_weight(e_prime, code_length) != n0_w)
    {
	printf("\nweight_of_e_decap \t %d\n",weight(e_prime, code_length));
        return -1;
    }

        

	m1 = (unsigned char *)malloc(BITS_TO_BYTES(k_prime));
   	 rho1 = (unsigned char *)malloc(BITS_TO_BYTES(k_sec));

    // Optimize modulo and removed copy to u1
    	memcpy(rho1, mot, BITS_TO_BYTES(k_sec));
    	memcpy(m1, mot + BITS_TO_BYTES(k_sec), BITS_TO_BYTES(code_dimension - k_sec));
    	free(mot);
	
	

        //Recover u1
        
	
	
	r1 = (unsigned char *)malloc(BITS_TO_BYTES(code_dimension));
	d1 = (unsigned char *)malloc(BITS_TO_BYTES(k_prime));
    // Compute r1 = G(m1) where G is composed of sponge SHA-512 function and extend function.
    // m_extend is no longer required because we are using KangarooTwelve which handles sizing

    // m: input type unsigned char len k_prime | r: output type unsigned char len code_dimesion
    test = KangarooTwelve(m1, BITS_TO_BYTES(k_prime), r1, BITS_TO_BYTES(code_dimension), custom, cus_len);
    assert(test == 0); // Catch Error

	

	// Compute d1 = H(m1) where H is  sponge SHA-512 function

    test = KangarooTwelve(m1, BITS_TO_BYTES(k_prime), d1, BITS_TO_BYTES(k_prime), custom, cus_len);
    assert(test == 0); // Catch Error

	

    // Return -1 if d distinct d1.
    // d starts at ct+code_length.
    if (memcmp(ct + BITS_TO_BYTES(code_length), d1, BITS_TO_BYTES(k_prime)) != 0)
    {
	printf("\ncompare_ct_d1\n");
        return -1;
    }
    free(d1);

	

       

     /*Step_5 of the decapasulation: Parse r1 as (rho2||sigma1).....................................................
      *********************************************************************************************************/


	 /*
   * Step_5 of the decapsulation: Parse r1 as (rho2||sigma1)
   */
    rho2 = (unsigned char *)malloc(BITS_TO_BYTES(k_sec));
    sigma = (unsigned char *)malloc(BITS_TO_BYTES(code_dimension));
    memcpy(rho2,r1,BITS_TO_BYTES(k_sec));
    memcpy(sigma,r1+BITS_TO_BYTES(k_sec),BITS_TO_BYTES(k_prime));

    //Return ⊥ if rho1 distinct rho2
    
    if (memcmp(rho1, rho2, BITS_TO_BYTES(k_sec)) != 0)
    {	
	printf("\ncompare_rho1_rho2\n");
        return -1;
    }
    free(r1);
    free(rho1);
    free(rho2);
	
    //Return ⊥ if rho1 distinct rho2
    



      
        //Generate error vector e2 of length code_length and weight n0_w from hash_sigma1 by using random_e function.


	hash_sigma = (unsigned char *)malloc(BITS_TO_BYTES(code_length));

    //Hashing sigma_extend by using KangarooTwelve function.

    test = KangarooTwelve(sigma, BITS_TO_BYTES(k_prime), hash_sigma, BITS_TO_BYTES(code_length), custom, cus_len);
    assert(test == 0); // Catch Error
    free(sigma);

    //Generate error vector e2 of length code_length and weight n0_w from
    //hash_sigma1 by using random_e function.
    e2 = random_e(code_length, 2, n0_w, hash_sigma);
    free(hash_sigma);



    /*
   * Step_7 of the decapsulation: Return ⊥ if e_prime distinct e.
   */
    if (memcmp(e_prime, e2, BITS_TO_BYTES(code_length)) != 0)
    {
	printf("\ncompare_e1_e2\n");
        return -1;
    }
    free(e_prime);
    free(e2);

    /*
   * Step_7 of the decapsulation: If the previous condition is not satisfied,
   * compute the shared secret ss by using KangarooTwelve 
   */
    test = KangarooTwelve(m1, BITS_TO_BYTES(k_prime), ss, ss_length, custom, cus_len);
    assert(test == 0); // Catch Error
    free(m1);

	printf("\n*************************Decapsulation:::La clef de session est:****************************\n\n");
            for(i = 0;i<ss_length;i++) printf("%d ",ss[i]);

            printf("\nFIN\n");		


	


    return 0;

}
/*END*/


