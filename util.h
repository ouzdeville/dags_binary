
#include<stdlib.h>
#include<stdio.h>
#include<sys/types.h>
#include<stdint.h>
#include"matrix.h"
#include"poly.h"
#include"time.h"



unsigned char* extend(unsigned char* m, int size1 , int size2);

int weight(unsigned char *r, int size);


unsigned char *random_m(int size);

int indice_in_vec(unsigned int * v, int j, int size);

unsigned char* random_e(int size, int q, int w,  unsigned char* sigma);


int  compare(unsigned char* tab1,unsigned char* tab2, int size);

unsigned char* gf_to_char(gf* a, int lenght);


gf sum_vect_element(gf* w,int length);
void set_vy_from_sk(gf* v, gf * y, const unsigned char * sk);
void set_sk_from_vy(gf * v, gf * y,unsigned char * sk);
void set_y_from_uvz(gf * u, gf* v, gf * Z, gf* y);
void set_vy_from_sk(gf* v, gf * y, const unsigned char * sk);
int bin_weight(unsigned char *r, int size);






