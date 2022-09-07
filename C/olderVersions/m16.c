#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

int msize = 32;
mpz_t *n, *nn, prec, tmp, eApprox;

void *m(long k, mpz_t* res){

    mpz_init_set_ui(n[0], k - 1);
    mpz_init_set_ui(n[msize-1], k + msize - 1);

    for(int i = 1; i < msize - 1; i++){
        mpz_set_ui(n[i], k + i);
    }
    for(int i = 1; i < msize-1; i+=2){
        int ind = (i+1) / 2;
        printf("ind:%d\n", ind);
        mpz_mul(nn[ind], n[i], n[i+1]);
    }

    mpz_set_ui(tmp, 1);
    for(int i = 1; i < msize / 2; i++){
        mpz_mul(tmp, tmp, nn[i]);
    }
    mpz_div(prec, prec, tmp);


    mpz_set_ui(tmp, 1);
    for(int i = 1, j = 2; i < msize / 2; i++, j+=2){
        mpz_mul(tmp, tmp, nn[i]);
        mpz_add(tmp, tmp, n[j]);
    }
    mpz_mul(*res, prec, tmp);

    mpz_add(eApprox, eApprox, *res);
    return NULL;
}

int main(){
    long accuracy = (int) pow(10., 4.);
    mpz_t res;
    mpz_inits(res, tmp, NULL);
    mpz_init_set_ui(prec, 10);
    mpz_pow_ui(prec, prec, (long) (accuracy + 2.*log2(accuracy)));
    mpz_set(eApprox, prec);
    
    n  = (mpz_t *) malloc(msize*sizeof(mpz_t));
    nn = (mpz_t *) malloc((int) (msize/2)*sizeof(mpz_t));
    for(int i = 0; i < msize; i++){
        mpz_init(n[i]);
    }
    for(int i = 1; i < msize / 2; i++){
        mpz_init(nn[i]);
    }


    clock_t t = clock();
    int i = 0;
    int k = 2;
    for(;;){
        m(k, &res);
        k += msize;
        i++;
        if(mpz_cmp_ui(prec, 0) == 0) break;
    }

    t = clock() - t;
    printf("Time taken %d:%d:%d (mm:ss:ms)\n", (int) t/60000, (int) (t/1000) % 60, t%1000);
    printf("Iterations: %d\n", i);
    
    
    FILE * fp = fopen("factorial.txt", "w+");
    mpz_out_str(fp, 10, eApprox);
    fclose(fp);

    system("\"C:/Program Files/nodejs/node\" ../Js/E/check.js ./factorial.txt ../e10M.txt");
}