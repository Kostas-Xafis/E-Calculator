#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <unistd.h>
#include <pthread.h>
// gcc calculateE7.c -lgmp -lpthread -o calc7
#define EXPO 7
#define E_BASE 10.05
#define PRECISION (int) pow(E_BASE, EXPO) // 10 million digits
#define BASE 10

int threads, done = 0, len_pow = 17;
long len, hlen, *nn_ind, split_len, split_hlen, **args;
mpz_t *n, *nn, eApprox, prec, tmp,
      *bfact_res, *fact, threads_res;
pthread_t tid;


void *bfact(long a, long b, long thread_ind) { //Binary splitting factorial
	if (!((a + b) % 2) && b - a > 2) {
        long m = (long) (a + b) / 2;
		mpz_t temp_res;

		bfact(a, m, thread_ind);
		mpz_init_set(temp_res, bfact_res[thread_ind]);
		bfact(m, b, thread_ind);
		mpz_mul(bfact_res[thread_ind], bfact_res[thread_ind], temp_res);

		mpz_clear(temp_res);
	} else {
		mpz_set_ui(fact[thread_ind], 1);
		for (long i = a + 1; i <= b; i++) mpz_mul_ui(fact[thread_ind], fact[thread_ind], i); 
		mpz_set(bfact_res[thread_ind], fact[thread_ind]); //bfact_res is also being reseted with no overlap here
        mpz_set(nn[nn_ind[thread_ind]], fact[thread_ind]);
        nn_ind[thread_ind]++;
	}
    return NULL;
}

void *start_bfact(void *vargv){
    long *args = (long *) vargv;
    long a = args[0];
    long b = args[1];
    long ind = args[2];
    bfact(a, b, ind);
    done++;
    return NULL;
}

void *m(long k, mpz_t *res){ 
    for(int i = 0; i < len; i++) mpz_set_d(n[i], k - 1 + i);

    for(int i = 0; i < threads; i++){
        args[i] = (long *) malloc(sizeof(long) * 3);
        args[i][0] = k+split_len*i-2;
        args[i][1] = k+split_len*(i+1)-2;
        args[i][2] = i;
        pthread_create(&tid, NULL, start_bfact, (void *)(args[i]));
    }
    pthread_join(tid, NULL);
    while(done != threads){}
    done = 0;
    for(int i = 0; i < threads; i++){
        mpz_mul(threads_res, threads_res, bfact_res[i]);
        nn_ind[i] = split_hlen*i;
    }
    

    mpz_div(prec, prec, threads_res);

    mpz_set_d(threads_res, 1);
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
    long i = 2;
    int times = 1;
    mpz_t rtmp;
    mpz_init(rtmp);

    time_t raw_time = time(NULL);
    printf("\nStart time: %s",ctime(&raw_time));
    printf("----Start----\n");
    
    while(1){
        mpz_set_d(threads_res, 1);
        m(i, &rtmp);
        mpz_add(eApprox, eApprox, rtmp);
        i+= len;
        if(i > 100000*times){
            printf("iterations: %ld\n", i);
            times++;
        }
        if(mpz_cmp_d(rtmp, 0) == 0) break;
    }
    mpz_clear(rtmp);
    return NULL;
}

int main(int argc, char const *argv[]){
        // Core variables
    threads = atoi(argv[argc-1]);
    len = (int) pow((double)2, (double)len_pow);
    hlen = len >> 1;
    split_len = len / threads; 
    split_hlen = split_len >> 1;

        // Allocating memory
    args = (long **) malloc(sizeof(long *)*threads);
    n = (mpz_t*) malloc(sizeof(mpz_t) * len);
    nn = (mpz_t*) malloc(sizeof(mpz_t) * hlen); //n(x)*n(x+1) array
    nn_ind = (long *) malloc(sizeof(long) * threads); //An array of indexes for threads
    bfact_res = (mpz_t*) malloc(sizeof(mpz_t)*threads); //Storing the results from the threads
    fact = (mpz_t*) malloc(sizeof(mpz_t)*threads);

        //Initializing variables
    for(int i = 0; i < len; i++) mpz_init(n[i]);
    for(int i = 0; i < hlen; i++) mpz_init(nn[i]);
    for(int i = 0; i < threads; i++){
        nn_ind[i] = split_hlen*i; //index offset of length of the threads split
        mpz_inits(bfact_res[i], fact[i], NULL);
    }
    mpz_inits(tmp, fact, NULL);
    mpz_init_set_d(threads_res, 1);
    mpz_init_set_d(prec, 10);
    mpz_pow_ui(prec, prec, PRECISION);
    mpz_init_set(eApprox, prec);


    clock_t t = clock();
    printf("----Initial Variables----\n length power: %d\n exponent base: %f\n exponent: %d\n approximate digits: %d\n threads: %d\n Split len: %ld\n Split hlen: %ld\n", len_pow, E_BASE, EXPO, PRECISION, threads, split_len, split_hlen);
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
    mpz_out_str(fp, 10, threads_res);
    fclose(fp);

        //Free memory
    free(n);
    free(nn);
    free(nn_ind);
    free(bfact_res);
    free(fact);
    free(args);
    mpz_clears(tmp, prec, eApprox, threads_res, NULL);
    pthread_exit(NULL);
        
    if(EXPO == 8)
        system("\"C:/Program Files/nodejs/node\" ../Js/E/check.js ./result.txt ../e100M.txt");
    else 
        system("\"C:/Program Files/nodejs/node\" ../Js/E/check.js ./result.txt ../e10M.txt");
    return 0;
}
