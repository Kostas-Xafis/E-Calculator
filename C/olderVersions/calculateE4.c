#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
// gcc -m32 calculateE4.c -lgmp -o calc4
#define PRECISION pow(10, 6) // 10 million digits
#define BASE 10

//Make a list u dumb ass

mpz_t eApprox, prec, tmp1,
      nm1, n1, n2, n3, n4, n5, n6, n7,
      n8, n9, n10, n11, n12, n13, n14, n15,
      n1n2, n3n4, n5n6, n7n8, n9n10, n11n12, n13n14;

void *nn(mpz_t* nn, mpz_t* a, mpz_t* b, long n){
    mpz_set_d(*a, n);
    mpz_set_d(*b, n+1);
    mpz_mul(*nn, *b, *a);
}

void *m16(long n, mpz_t* res){
    mpz_set_ui(nm1, n - 1);
    mpz_set_ui(n15, n+15);    
	nn(&n1n2, &n1, &n2, n+1);
    nn(&n3n4, &n3, &n4, n+3);
	nn(&n5n6, &n5, &n6, n+5);
	nn(&n7n8, &n7, &n8, n+7);
	nn(&n9n10, &n9, &n10, n+9);
	nn(&n11n12, &n11, &n12, n+11);
    nn(&n13n14, &n13, &n14, n+13);
	
	mpz_mul_ui(tmp1, nm1, n);
    mpz_mul(tmp1, tmp1, n1n2);
    mpz_mul(tmp1, tmp1, n3n4);
    mpz_mul(tmp1, tmp1, n5n6);
    mpz_mul(tmp1, tmp1, n7n8);
    mpz_mul(tmp1, tmp1, n9n10);
    mpz_mul(tmp1, tmp1, n11n12);
    mpz_mul(tmp1, tmp1, n13n14);
    mpz_div(prec, prec, tmp1);

    mpz_mul(tmp1, n1n2, n1);
    mpz_add(tmp1, tmp1, n3);
    mpz_mul(tmp1, tmp1, n3n4);
    mpz_add(tmp1, tmp1, n5);
    mpz_mul(tmp1, tmp1, n5n6);
    mpz_add(tmp1, tmp1, n7);
    mpz_mul(tmp1, tmp1, n7n8);
    mpz_add(tmp1, tmp1, n9);
    mpz_mul(tmp1, tmp1, n9n10);
    mpz_add(tmp1, tmp1, n11);
    mpz_mul(tmp1, tmp1, n11n12);
    mpz_add(tmp1, tmp1, n13);
    mpz_mul(tmp1, tmp1, n13n14);
    mpz_add(tmp1, tmp1, n15);
    mpz_mul(*res, prec, tmp1);

    return NULL;
}

void *calculation(){

    mpz_t tmp;
    long i = 2;
    `mpz_init(tmp);
    printf("----Start----\n");
    while(1){
        m16(i, &tmp);
        mpz_add(eApprox, eApprox, tmp);
        i+= 16;
        if(mpz_cmp_d(tmp, 0) == 0) break;
    }
    return NULL;
}

int main(int argc, char const *argv[]){
        //Initializing variables
    mpz_init(nm1, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, 
        n11, n12, n13, n14, n15, n1n2, n3n4, n5n6, n7n8, n9n10,
        n11n12, n13n14, tmp1, NULL);
    mpz_init_set_d(prec, 10);
    mpz_pow_ui(prec, prec, PRECISION);
    mpz_init_set(eApprox, prec);

    clock_t start = clock(), diff;
    calculation();
        //Perfomance
    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("----Finished----\n");
    printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);

    char* str = mpz_get_str(NULL, BASE, eApprox);
    mpz_clear(eApprox);
    mpz_clear(prec);

    FILE *fp;
    fp = fopen("result.txt", "w");
    fwrite(str, strlen(str), sizeof(char), fp);
    fclose(fp);

    system("node ../Js/E/check.js ./result.txt ../e.txt");
    return 0;
}
