#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <unistd.h>
#include <pthread.h>
// #include <signal.h>
// gcc -m32 calculateE3.c -lgmp -lpthread -o calc3
#define PRECISION pow(10, 6) // 10 million digits
#define HPRECISION PRECISION/2
#define BASE 10
typedef struct{
    int expo2; // Exponent of 2
    mpz_t rem; // Remainder
}factors;

factors* factorArray;
mpz_t finalApprox;
clock_t endTime, start;
int done = 0,
    calculating = 0,
    threads,
    pause = 0;

void findFactors(){
    mpz_t rem;
    mpz_init(rem);
    for(int i = 0; i < HPRECISION; i++){
        mpz_set_d(rem, i);
        for(int k = i+1; k < i+threads; k++){mpz_mul_ui(rem, rem, k);}
        int expo2 = 0;
        while(mpz_divisible_ui_p(rem, 2) && mpz_cmp_ui(rem, 0) ){
            mpz_div_2exp(rem, rem, 1);
            expo2++;
        }
        factorArray[i].expo2 = expo2;
        mpz_init_set(factorArray[i].rem, rem);
    }
    mpz_init_set_d(factorArray[0].rem, 1);
}


void *calculation(void *vargv){
        //Initializing variables
    mpz_t eApprox, prec;
    mpz_init_set_d(prec, 10);
    mpz_pow_ui(prec, prec, PRECISION);
    mpz_init_set_ui(eApprox, 0);

    int iter = PRECISION/2,
        times = 1,
        i = (int)(vargv) + 1, //offset + 1
        expo2 = 1;
    for(int k=1; k < i; expo2*=k++){};
    mpz_fdiv_q_ui(prec, prec, expo2);
    
    for (; i <= iter; i+= threads){
        while(pause == 1){ //Pausing the execution
            sleep(1);
        }
        mpz_add(eApprox, eApprox, prec);
        
        expo2 = factorArray[i].expo2;
        if(expo2)
            mpz_div_2exp(prec, prec, expo2);
        mpz_fdiv_q(prec, prec, factorArray[i].rem);

        if(mpz_cmp_ui(prec, 0) == 0) break;
        // if (i >= 5000 * times) {
		// 	printf("Iterations: %d\n", 5000 * times);
		// 	times++;
		// }
    }
        //prevent from over lapping
    while(calculating == 1) sleep(0.01);
    calculating = 1;
    mpz_add(finalApprox, finalApprox, eApprox);
    calculating = 0;
    
    done++;
    if(done == threads)
        endTime = clock();
    
    mpz_clear(eApprox);
    mpz_clear(prec);
    return NULL;
}

void *PauseCalc(){
        //Perfomance
    int paused = 0, input;
    clock_t time = 0, pauseTime;
    while(done != threads){
        // paused = pause;
        sleep(1);
        // scanf("%d", &input);
        // if(input == 0 || input == 1)//Pause logic
        //     pause = input;
        // if(pause == 1 && paused == 0)//Performance logic
        //     pauseTime = clock();
        // if(pause == 0 && paused == 1)
        //     time += clock() - pauseTime;
        // //!Need to remove the time from the last input;
    }
    time = endTime - start - time;
    int msec = time * 1000 / CLOCKS_PER_SEC;
    printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
}

    // Arguments #1: Number of threads
int main(int argc, char const *argv[]){
    threads = atoi(argv[argc-1]);
    factorArray = (factors *) malloc(sizeof(factors) * HPRECISION);
    findFactors();
    
    mpz_t p;
    mpz_init_set_d(p, 10);
    mpz_pow_ui(p, p, PRECISION);
    mpz_init_set(finalApprox, p);
    mpz_clear(p);

        //Create Threads
    pthread_t tid;
    for(int i = 0; i < threads; i++){
        pthread_create(&tid, NULL, calculation, (void *)i);
    }
    pthread_create(&tid, NULL, PauseCalc, "");
    printf("Threads Created\n");
    start = clock();
    pthread_join(tid, NULL);

        //Store to file
    FILE *fp;
    fp = fopen("result.txt", "w");
    mpz_out_str(fp, BASE, finalApprox);
    fclose(fp);

    mpz_clear(finalApprox);
    free(factorArray);
        //Compare the result
    system("node ../Js/E/check.js ./result.txt ../e.txt");

    return 0;
}
