#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <unistd.h>
#include <pthread.h>
// gcc -m32 calculateE2.c -lgmp -lpthread -o calc2
#define PRECISION pow(10, 6) // 10 million digits
#define BASE 10

mpz_t finalApprox;
int done = 0,
    calculating = 0,
    threads;
void *calculation(void *vargv){
        //Initializing variables
    mpz_t eApprox, prec, sub_fact_increment, sub_fact_mul;
    mpz_init_set_d(prec, 10);
    mpz_pow_ui(prec, prec, PRECISION);
    mpz_init_set_d(sub_fact_increment, 1);
    mpz_init_set_d(eApprox, 0);
    mpz_init_set_d(sub_fact_mul, 1);

    int iter = PRECISION/2,
        times = 1, 
        k = 1, 
        i = (int)(vargv) + 1; //offset + 1

    for(; k < i; k++){
        mpz_add_ui(sub_fact_increment, sub_fact_increment, 1);
        mpz_mul(sub_fact_mul, sub_fact_mul, sub_fact_increment);
    }

    for (; i <= iter; i+= threads){
        mpz_fdiv_q(prec, prec, sub_fact_mul);
        mpz_add(eApprox, eApprox, prec);

        mpz_set_d(sub_fact_mul, 1);
        for(k = 1; k <= threads; k++){
            mpz_add_ui(sub_fact_increment, sub_fact_increment, 1);
            mpz_mul(sub_fact_mul, sub_fact_mul, sub_fact_increment);
        }

        if(mpz_cmp_ui(prec, 0) == 0) break;
        if (i >= 5000 * times) {
			printf("Iterations: %d\n", 5000 * times);
			times++;
		}
    }
        //prevent from over lapping
    while(calculating == 1) sleep(0.01);
    calculating = 1;
    mpz_add(finalApprox, finalApprox, eApprox);
    calculating = 0;

    mpz_clear(eApprox);
    mpz_clear(prec);
    mpz_clear(sub_fact_increment);
    mpz_clear(sub_fact_mul);
    done++;
    return NULL;
}

int main(int argc, char const *argv[]){
    threads = atoi(argv[argc-1]);
    mpz_t p;
    mpz_init_set_d(p, 10);
    mpz_pow_ui(p, p, PRECISION);
    mpz_init_set(finalApprox, p);
    mpz_clear(p);

    pthread_t tid;
    for(int i = 0; i < threads; i++){
        pthread_create(&tid, NULL, calculation, (void *)i);
    }
    printf("Threads Created\n");
    clock_t start = clock(), diff;
    pthread_join(tid, NULL);
        //Perfomance
    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);

    char* str = mpz_get_str(NULL, BASE, finalApprox);

    FILE *fp;
    fp = fopen("result.txt", "w");
    fwrite(str, strlen(str), sizeof(char), fp);
    fclose(fp);

    mpz_clear(finalApprox);

    system("node ../Js/check.js");

    return 0;
}
