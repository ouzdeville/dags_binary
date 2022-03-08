/******************************************************************************************
* BINARY DAGS: Key Encapsulation using binary Dyadic GS Codes.                            *
* This is the binary version of DAGS  submitted to the NIST Post=Quantum Cryptography.    *
* For any other usage , contact the author(s) to ascode_dimmension permission.            *
*******************************************************************************************
*/

#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <stdint.h>
#include <assert.h>
#include <keccak/KangarooTwelve.h>

#include "round.h"

//#include "param.h"

#include "util.h"

#define cus_len 4 // I random pick this number to fullfill parameter

/**
 * @brief Method to perform the encapsulation using public key
 * @param pk public key from the scheme
 * @param ct
 * @param ss secret shared
 */
int encapsulation(const unsigned char *pk, unsigned char *ct, unsigned char *ss);
