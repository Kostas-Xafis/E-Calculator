#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
// gcc -m32 calculateE.c -lgmp -o calc
#define PRECISION pow(10, 6) // 10 million digits
#define HPRECISION PRECISION/2
#define BASE 10
typedef struct{
    int exp2; // Exponent of 2
    int rem; // Remainder
}factors;

factors* factorArray;
void *findFactors(){
    for(int i = 0; i < HPRECISION; i++){
        int rem = i;
        int exp2 = 0;
        // int digits = (int) floor(log2(i));
        while(!(rem % 2) && rem != 0){
            rem = rem >> 1;
            exp2++;
        }
        factorArray[i].exp2 = exp2;
        factorArray[i].rem = rem;
    }
    factorArray[0].rem = 1;
}

void *calculation(mpz_t* eApprox, mpz_t* prec){
    int iter = HPRECISION, times = 1;
    int exp2 = 0;
    for (int i = 0; i < iter; i++){
        exp2 = factorArray[i].exp2;
        if(exp2)
            mpz_tdiv_q_2exp(*prec, *prec, exp2);
        mpz_fdiv_q_ui(*prec, *prec, factorArray[i].rem);
        // mpz_fdiv_q(*prec, *prec, *sub_fact);
        mpz_add(*eApprox, *eApprox, *prec);
        // mpz_add_ui(*sub_fact, *sub_fact, 1);
        if(mpz_cmp_ui(*prec, 0) == 0) break;
        if (i >= 5000 * times) {
			printf("Iterations: %d\n", 5000 * times);
			times++;
		}
    }
    return NULL;
}

int main(int argc, char const *argv[]){
    factorArray = (factors *) malloc(sizeof(factors) * HPRECISION);
    findFactors();
        //Initializing variables
    mpz_t eApprox, prec;
    mpz_init_set_d(prec, 10);
    mpz_pow_ui(prec, prec, PRECISION);
    mpz_init_set_d(eApprox, 0);

    clock_t start = clock(), diff;
    calculation(&eApprox, &prec);
        //Perfomance
    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);

    char* str = mpz_get_str(NULL, BASE, eApprox);
    mpz_clear(eApprox);
    mpz_clear(prec);
    free(factorArray);

    FILE *fp;
    fp = fopen("result.txt", "w");
    fwrite(str, strlen(str), sizeof(char), fp);
    fclose(fp);

    system("node ../Js/check.js");
    return 0;
}
