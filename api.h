
#ifndef __API_H_INCLUDED__
#define __API_H_INCLUDED__

#define CRYPTO_SECRETKEYBYTES 2*code_length*sizeof(gf) //4*code_length //4418 //   8642  //3313.5
#define CRYPTO_PUBLICKEYBYTES   BITS_TO_BYTES((code_dimension*mst))//(BITS_TO_LONG(code_length-code_dimension) * sizeof(long) * code_dimension)
#define CRYPTO_CIPHERTEXTBYTES (BITS_TO_BYTES(code_length)+BITS_TO_BYTES(k_prime))   // 1616
#define CRYPTO_BYTES 64

#define CRYPTO_ALGNAME "DAGS_128"



int crypto_kem_keypair(
    unsigned char *pk,
    unsigned char *sk);

int crypto_kem_enc(
    unsigned char *ct,
    unsigned char *ss,
    const unsigned char *pk);

int crypto_kem_dec(
    unsigned char *ss,
    const unsigned char *ct,
    const unsigned char *sk);





#endif
