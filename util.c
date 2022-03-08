/******************************************************************************************
* BINARY DAGS: Key Encapsulation using binary Dyadic GS Codes.                            *
* This is the binary version of DAGS  submitted to the NIST Post=Quantum Cryptography.    *
* For any other usage , contact the author(s) to ascode_dimmension permission.            *
*******************************************************************************************
*/


#include"util.h"
#include"gf.h"
#include"matrix.h"
#include <string.h>
#include"api.h"




unsigned char* extend(unsigned char* m, int size1 , int size2){
   int i;
            unsigned char* res = (unsigned char*)calloc(size2,sizeof(unsigned char));
            for (i=0;i< size1;i++) res[i] = m[i];
return res;
}


/*weight computes the weight of a sequence of elements of type unsigned char***********************
  *************************************************************************************************/


int weight(unsigned char *r, int size)
{
	int i = 0, w = 0;
	for (i = 0; i < size; i++)
	{
		if (r[i] != 0)
			w++;
	}
	return w;
}

int bin_weight(unsigned char *r, int size)
{
	int i = 0, w = 0;
	for (i = 0; i < size; i++)
	{
		if ((r[i/8]>>(i%8))&1)
			w++;
	}
	return w;
}

//*END**********************************************************************************************/


/*random_m generate randomly a sequence of size element of F_q*************************************
  *************************************************************************************************/

unsigned char *random_m(int size)
{
	unsigned char *r = (unsigned char *)malloc(size);
	int i;
	randombytes(r, size);
	return r;
}


/*indice_in_vec test if element is in tab **********************************************************
  *************************************************************************************************/

int indice_in_vec(unsigned int * v, int j, int size){
   int i;
for(i=0;i<size;i++){
            if(v[i] == j) return 1;

  }
return 0;

}
//*END**********************************************************************************************/

/*random_e **************************************************************************************
  *************************************************************************************************/

unsigned char* random_e(int size, int q, int w,  unsigned char* sigma)
{
 
	unsigned char *e = (unsigned char *)calloc(BITS_TO_BYTES(size), sizeof(unsigned char));
	unsigned int *v = (unsigned int *)calloc(size, sizeof(unsigned int));
	int j = 0, k = 0, jeton = 0;

	while (1)
	{


		do
		{
			jeton = (sigma[k + 1] ^ (sigma[k] << 4)) % size;

			k++;
		} while (indice_in_vec(v, jeton, j + 1) == 1); //Only check j elements
		v[j] = jeton;
		e[jeton/8] |= (1ULL << jeton%8);
		jeton = 0;
		j++;
		if (j == w)
				{
					break;
				}
	}
	free(v);
	return e;
}
//*END**********************************************************************************************/



/*compare******************************************************************************************
 *************************************************************************************************/
int  compare(unsigned char* tab1,unsigned char* tab2, int size){
         int i=0;
         for(i=0;i<size;i++){
               if(tab1[i] != tab2[i]) return 0;
         }
         return 1;
}





//*END**********************************************************************************************/

/**
 *
 */
gf sum_vect_element(gf* w,int length){
   gf tmp=0;
   int j=0;
        for (j = 0; j < length; j++)
        {
            tmp ^= w[j];
        }
   return tmp;
}


void set_sk_from_vy(gf * v, gf * y,unsigned char * sk)
{
	memcpy(sk, v, code_length*sizeof(gf_t));
	memcpy(sk+code_length*sizeof(gf_t), y, sizeof(gf_t)*code_length);

}

void set_y_from_uvz(gf * u, gf* v, gf * Z, gf* y){

        gf* z = (gf *)calloc(code_length, sizeof(gf));
        gf pol,aux;
	int i,j;

	for (i = 0; i < n0_val; i++)
	    {
		for (j = 0; j < order; j++)
		{
		    z[i * (order) + j] = Z[i];
		}
	    }
	
        for (i = 0; i <  code_length; i++)
        {
            pol =1;
            for (j = 0; j < order; j++)
            {
		aux = v[i]^u[j];
                pol = gf_mul(pol,aux);
            }
	y[i] = gf_mul(z[i], gf_pow(gf_inv(pol),pol_deg));

        }
        free(z);
}



void set_vy_from_sk(gf* v, gf * y, const unsigned char * sk)
{
	memcpy(v, sk, code_length*sizeof(gf_t));
		//sk += code_length*sizeof(gf_t);
	memcpy(y,sk+code_length*sizeof(gf_t), sizeof(gf_t)*code_length);
	
	
}






