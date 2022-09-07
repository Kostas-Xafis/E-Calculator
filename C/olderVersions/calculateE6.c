#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
// gcc -m32 calculateE6.c -lgmp -o calc6
#define EXPO 7
#define E_BASE 11
#define PRECISION (int) pow(E_BASE, EXPO) // 10 million digits
#define BASE 10


//make a hashtable for the multiplication of the leafs of bfact
long len, hlen;
mpz_t *n, *nn, eApprox, prec, tmp,
       bfact_res, fact;
int nn_ind = 0;

void* bfact(long a, long b) { //Binary splitting factorial
	if (!((a + b) % 2) && b - a > 2) {
        long m = (long) (a + b) / 2;
		mpz_t temp_res;

		bfact(a, m);
		mpz_init_set(temp_res, bfact_res);
		bfact(m, b);
		mpz_mul(bfact_res, bfact_res, temp_res);

		mpz_clear(temp_res);
	} else {
		mpz_set_ui(fact, 1);
		for (long i = a + 1; i <= b; i++) mpz_mul_ui(fact, fact, i); 
		mpz_set(bfact_res, fact); //bfact_res is also being reseted safely here
        mpz_set(nn[nn_ind], fact);
        nn_ind++;
	}
}

void *m(long k, mpz_t* res){
    for(int i = 0; i < len; i++)
        mpz_set_d(n[i], k - 1 + i);

    bfact(k-2, k + len - 2);
    mpz_div(prec, prec, bfact_res);
    nn_ind = 0;

    mpz_set_d(tmp, 0);
    for(int i = 1; i < hlen; i++){
        mpz_add(tmp, tmp, n[2*i]);
        mpz_mul(tmp, tmp, nn[i]);
    }
    mpz_add_ui(tmp, tmp, k+len-1);
    mpz_mul(*res, prec, tmp);

    return NULL;
}

void *calculation(){
    mpz_t rtmp;
    long i = 2;
    int times = 1;
    mpz_init(rtmp);

    time_t raw_time = time(NULL);
    printf("\nStart time: %s",ctime(&raw_time));
    printf("----Start----\n");
    
    while(1){
        m(i, &rtmp);
        // if(i >= 5*len) abort();
        mpz_add(eApprox, eApprox, rtmp);
        i+= len;
        if(i > 100000*times){
            printf("iterations: %d\n", i);
            times++;
        }        
        if(mpz_cmp_d(rtmp, 0) == 0) break;
    }
    mpz_clear(rtmp);
    return NULL;
}

int main(int argc, char const *argv[]){
    int len_pow = 16;
    len = (int) pow((double)2, (double)len_pow);
    hlen = len >> 1;
    n = (mpz_t*) malloc(sizeof(mpz_t) * len);
    nn = (mpz_t*) malloc(sizeof(mpz_t) * hlen);
        //Initializing variables
    for(int i = 0; i < len; i++) mpz_init(n[i]);
    for(int i = 0; i < hlen; i++) mpz_init(nn[i]);
    mpz_inits(tmp, bfact_res, fact, NULL);
    mpz_init_set_d(prec, 10);
    mpz_pow_ui(prec, prec, PRECISION);
    mpz_init_set(eApprox, prec);

    clock_t t = clock();
    printf("----Initial Variables----\n length power: %d\n exponent base: %d\n exponent: %d\n approximate digits: %d\n", len_pow, E_BASE, EXPO, PRECISION);
    calculation();
        
        //Perfomance
    t = clock() - t;
    int mills = t * 1000 / CLOCKS_PER_SEC;
    int sec = mills/1000;
    int min = sec / 60;
    printf("----Finished----\n");
    printf("Time taken %d:%d:%d (mm:ss:ms)\n", min, sec%60, mills%1000);
    
        //Store
    FILE *fp;
    fp = fopen("result.txt", "w");
    mpz_out_str(fp, 10, eApprox);
    fclose(fp);

        //Free memory
    free(n);
    free(nn);
    mpz_clears(tmp, prec, eApprox, NULL);
        
        //Verify digits
    system("node ../Js/E/check.js ./result.txt ../e.txt");
    return 0;
}
