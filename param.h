#ifndef PARAM_H
#define PARAM_H
#define BITS_PER_LONG (8 * sizeof (unsigned long))
#define BITS_TO_BYTES(nb_bits) (((nb_bits) - 1) / 8 + 1)
#define pol_deg 2
#define code_length 6400
#define order 128
#define EXT_DEGREE 13
#define mst (order*pol_deg*EXT_DEGREE)
#define code_dimension (code_length-mst)
#define k_prime  256
#define k_sec (code_dimension-k_prime)
#define n0_val (code_length/order)
#define n0_w (order*pol_deg)/2
#define ss_length 64
#define debug 0
#endif
