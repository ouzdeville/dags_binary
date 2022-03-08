/******************************************************************************************
* BINARY DAGS: Key Encapsulation using binary Dyadic GS Codes.                            *
* This is the binary version of DAGS  submitted to the NIST Post=Quantum Cryptography.    *
* For any other usage , contact the author(s) to ascode_dimmension permission.            *
*******************************************************************************************
*/
#include "decoding.h"

//Bulding of decoding fuction

/*
 * The polynome_syndrome_1 function compute the syndrome S in polynomial form
 *  with inputs a short IV and a Parity matrix H
 */



void polynome_syndrome_1(binmat_t H, const unsigned char *mot, poly_t S)
{
    int i, j;

    gf tmp;
    for (j = 0; j < H.rown; j++)
    {
        tmp = 0;
        for (i = 0; i < H.coln; i++)
        {
            tmp ^= gf_mul(mat_coeff(H, j, i), ((gf)mot[i]));
        }
        S->coeff[j] = tmp;
    }
    poly_calcule_deg(S);
}



/*
 * The polynome_syndrome_from_vy function compute the syndrome S in polynomial form
 *  with inputs a short IV and two secret vector
 */
void polynome_syndrome_from_vy(gf* v, gf* y, const unsigned char *mot, poly_t S)
{
    int i, j;
    int st = (order)*pol_deg;
	

    gf* w = (gf *)calloc(code_length, sizeof(gf));
    for (j = 0; j < code_length; j++)
    {
         w[j]=gf_mul(y[j], ((gf)mot[j/8]>>(j%8))&1);
    }
    S->coeff[0]=sum_vect_element(w,code_length);
    
     for (i = 1; i < st; i++)
        {
		for (j = 0; j < code_length; j++)
		    {
			 w[j]=gf_mul(v[j], w[j]);
		    }	

	S->coeff[i]=sum_vect_element(w,code_length);            
        }
    poly_calcule_deg(S);
	
free(w);
}



int decoding_from_vy(gf* v,gf* y, const unsigned char *c, unsigned char *error,
               unsigned char *code_word)
{
    
    gf_init(EXT_DEGREE);
    int i, dr;
  
    int st = (order) * pol_deg;
    poly_t Syndrome;
    poly_t omega, sigma, re, uu, u, quotient, resto, app, temp;
    poly_t pol;
    
   

    //Compute Syndrome normally
    Syndrome = poly_alloc(st - 1);
    polynome_syndrome_from_vy(v,y, c, Syndrome);

    

    if (Syndrome->deg == -1)
    {

        return -1;
    }

    //Resolution of the key equation
    re = poly_alloc(st);
    re->coeff[st] = 1;
    poly_calcule_deg(re);
    resto = poly_alloc(st);
    resto->coeff[st] = 1;
    poly_calcule_deg(resto);

    app = poly_alloc(st);
    uu = poly_alloc(st);
    u = poly_alloc(st);
    //poly_set_to_unit(u);
    //Set to unit
    u->coeff[0] = 1;
    u->deg = 0;

    dr = Syndrome->deg;


    while (dr >= (st / 2))
    {
        quotient = poly_quo(re, Syndrome);
        poly_rem(resto, Syndrome);

        poly_set(re, Syndrome);
        poly_set(Syndrome, resto);
        poly_set(resto, re);

        poly_set(app, uu);
        poly_set(uu, u);

        temp = poly_mul(u, quotient);
        poly_free(u);
        poly_free(quotient);
        u = temp;
        poly_add_free(u, u, app);
        poly_calcule_deg(Syndrome);
        dr = Syndrome->deg;
    }

    poly_free(re);
    poly_free(uu);
    poly_free(app);
    poly_free(resto);

    //Then we find error locator poly (sigma) and error evaluator poly (omega)
    pol = poly_alloc(0);
    pol->coeff[0] = gf_inv(poly_eval(u, 0));
    omega = poly_mul(Syndrome, pol);
    poly_free(Syndrome);
    sigma = poly_mul(u, pol);
    poly_free(pol);
    poly_free(u);
  
	
	
	int d=poly_calcule_deg(sigma);
	
	
	
	
	// BTA: roots search algorithm
	gf_t res[n0_w];
	int rt = BTA(sigma, res);
	
	// Tab_ind: Array containing the position of each element into the support v
	int *Tab_ind;
	Tab_ind = (int *) calloc(gf_card(), sizeof(int));
	int pow=0;
	for (i = 0; i < code_length; i++){
		Tab_ind[gf_log[v[i]]]= i+1;
	}
	// Computing the error vector 
	int pos_error=0;
	for (i = 0; i < n0_w; i++){
		pow = gf_log[gf_inv(res[i])];
		pos_error = Tab_ind[pow]-1 ; 
		error[pos_error/8]|=1<<(pos_error%8);
		
	}
	

    //Reconstruction of code_word
    for (i = 0; i < code_length/8; i++)
    {
        
	code_word[i] = (c[i] ^ error[i]);
        
    }

	
     
    return 1;
}


