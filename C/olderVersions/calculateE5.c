#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
// gcc -m32 calculateE5.c -lgmp -o calc5
#define PRECISION (int) pow(10.15, 7) // 10 million digits 
#define BASE 10
//I think i must continue with lower length n-array to accurately measure 10^k
// instead of having that 10.15

long len, hlen;
mpz_t *n, *nn, eApprox, prec, tmp;

void *nn_calc(mpz_t* nn, mpz_t* a, mpz_t* b, long n){
    mpz_set_d(*a, n);
    mpz_set_d(*b, n+1);
    mpz_mul(*nn, *b, *a);
}

void *m16(long k, mpz_t* res){    
    for(int i = 0; i < hlen; i++)
        nn_calc(&nn[i], &n[2*i], &n[2*i+1], k + 2*(i-1) + 1);

        //Next up binary split
    mpz_set_d(tmp, 1);
    for(int i = 0; i < hlen; i++) mpz_mul(tmp, tmp, nn[i]);
    mpz_div(prec, prec, tmp);

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
    printf("\n----Start----\n");
    
    while(1){
        m16(i, &rtmp);
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
    len = (int) pow((double)2, (double)15);
    hlen = (int) len/2;
    n = (mpz_t*) malloc(sizeof(mpz_t) * len);
    nn = (mpz_t*) malloc(sizeof(mpz_t) * hlen);
        //Initializing variables
    for(int i = 0; i < len; i++) mpz_init(n[i]);
    for(int i = 0; i < hlen; i++) mpz_init(nn[i]);
    mpz_init(tmp);
    mpz_init_set_d(prec, 10);
    mpz_pow_ui(prec, prec, PRECISION);
    mpz_init_set(eApprox, prec);

    clock_t t = clock();
    calculation();
        
        //Perfomance
    t = clock() - t;
    int mills = t * 1000 / CLOCKS_PER_SEC;
    int sec = mills/1000;
    int min = sec / 60;
    printf("----Finished----\n\n");
    printf("Time taken %d:%d:%d (mm:ss:ms)\n", min, sec%60, mills%1000);
    
        //Store
    FILE *fp;
    fp = fopen("result.txt", "w");
    mpz_out_str(fp, 10, eApprox);
    fclose(fp);

        //Free memory
    free(n);
    free(nn);
    mpz_clear(tmp);
    mpz_clear(prec);
    mpz_clear(eApprox);

    system("node ../Js/E/check.js ./result.txt ../e.txt");
    return 0;
}
