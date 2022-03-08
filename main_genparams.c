
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
/*
13 6400 128 2
14 11520 256 2
14 14080 256 2
m   n    s   t
4   32   1   2
s=
mst < n
 */
int main(int argc, char ** argv) {
	int m, n, s, t;
	FILE * fichier;
	double * res, *res2, wf;

	m = n = s = t = 0;
	if (argc > 4) {
		m = atoi(argv[1]);
		n = atoi(argv[2]);
		s = atoi(argv[3]);
		t = atoi(argv[4]);
	}


	if ((m <= 0) || (n <= 0) || (s <= 0) || (t <= 0)) {
		fprintf(stderr, "Usage: %s m n s t [extension length order numBocks]\n",
				argv[0]);
		fprintf(stderr,
				"all arguments are positive integers, with m > 5, and 0 < t < 2^m/m\n");
		fprintf(stderr,
				"Look at the documentation for more information on the arguments\n");
		exit(0);
	}
	if (m*s*t>=n) {
				fprintf(stderr, "The parameters s, t for the GS code should be chosen such that mst < n.\n");
				exit(0);
			}
	if (n%s!=0) {
					fprintf(stderr, "The parameter s should divide n and be power of 2\n");
					exit(0);
	}

	fichier = fopen("param.h", "w");
	fprintf(fichier, "#ifndef PARAM_H\n");
	fprintf(fichier, "#define PARAM_H\n");

	fprintf(fichier, "#define BITS_PER_LONG (8 * sizeof (unsigned long))\n");
	fprintf(fichier, "#define BITS_TO_BYTES(nb_bits) (((nb_bits) - 1) / 8 + 1)\n");

	fprintf(fichier, "#define pol_deg %d\n", t);
	fprintf(fichier, "#define code_length %d\n", n);
	fprintf(fichier, "#define order %d\n", s);
	fprintf(fichier, "#define EXT_DEGREE %d\n", m);
	fprintf(fichier, "#define mst (order*pol_deg*EXT_DEGREE)\n");
	fprintf(fichier, "#define code_dimension (code_length-mst)\n");
	fprintf(fichier, "#define k_prime  256\n");
	fprintf(fichier, "#define k_sec (code_dimension-k_prime)\n");
	fprintf(fichier, "#define n0_val (code_length/order)\n");
	fprintf(fichier, "#define n0_w (order*pol_deg)/2\n");
	fprintf(fichier, "#define ss_length 64\n");
	fprintf(fichier, "#define debug 0\n");

	fprintf(fichier, "#endif\n");
	fclose(fichier);

	return 0;
}
