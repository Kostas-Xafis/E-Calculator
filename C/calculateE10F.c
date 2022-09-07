#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <unistd.h>
#include <pthread.h>

/*
gcc calculateE10F.c -I/c/msys64/gmp -L/c/msys64/gmp/lib -lgmp -lpthread -o calc10F
*/
#define EXPO 5
#define BASE 10.05
#define PRECISION (int) pow(BASE, EXPO)
#define THREADS 8
#define LEN_POW 15

int bfact_done = 0, final_thread_mul = 0, done_calc = 0, *threads_done, 
    queueSize = 20, qfront = 0, qrear = 0, fullQueue = 0, ResultPrecison, FactorialPrecision;
long len, hlen, *nn_ind, split_len, split_hlen, **args;
mpf_t *n, *nn, eApprox, *tmp,
      *bfact_res, *fact, *tmpQueue, *bfrQueue;
pthread_t tid;
clock_t t;

int mpf_t_size = sizeof(mpf_t);
void *pointerSwap(mpf_t *a, mpf_t *b){
    mpf_t *temp = (mpf_t *) malloc(mpf_t_size); 

    memcpy(temp, b, mpf_t_size);
    memcpy(b, a, mpf_t_size);
    memcpy(a, temp, mpf_t_size);
}

void *bfact(long a, long b, long thread_ind) { //Binary splitting factorial
	if (b - a > 2) {
        long m = (long) floor((a + b) / 2);
		mpf_t temp_res;
        mpf_init2(temp_res, FactorialPrecision);
		
        bfact(a, m, thread_ind);
		mpf_set(temp_res, bfact_res[thread_ind]);
		bfact(m, b, thread_ind);
		mpf_mul(bfact_res[thread_ind], bfact_res[thread_ind], temp_res);

		mpf_clear(temp_res);
	} else {
		mpf_set_ui(fact[thread_ind], 1);
		for (long i = a + 1; i <= b; i++) mpf_mul_ui(fact[thread_ind], fact[thread_ind], i); 
		mpf_set(bfact_res[thread_ind], fact[thread_ind]); //bfact_res is also being reseted with no overlap here
        mpf_set(nn[nn_ind[thread_ind]], fact[thread_ind]);
        nn_ind[thread_ind]++;
	}
    return NULL;
}

void *start_bfact(void * vargv){
    long *args = (long *) vargv;
    long a = args[0];
    long b = args[1];
    long ind = args[2];
    long k = args[3];

    bfact(a, b, ind);

    int sibling = 1;
    if(ind != 0) bfact_done++;
    if ((ind+1) % 2 == 0) {
        //Each odd indexed thread multiplies it's values with his left sibling's value
        //recursively until there is only 1 thread remaining.
        //ie:(t0*t1),(t2*t3),(t4*t5),(t6*t7) -> (t1*t3), (t5*t7) -> (t3*t7) -> return t7
        while (1) {
            if (threads_done[ind - sibling]) { //if left sibling thread is done
                mpf_mul(bfact_res[ind], bfact_res[ind], bfact_res[ind-sibling]);
                sibling *= 2;
                threads_done[ind] = 1; //prevent overlapping
                if (((ind+1) / sibling) % 2 == 1) break;
                threads_done[ind] = 0;
            }
        }
    } else threads_done[ind] = 1;

    if (ind == 0) {
        //Make the first thread
        for(int i = 0; i < len; i++) mpf_set_ui(n[i], k - 1 + i);
        while(bfact_done != THREADS-1){}
        for (int i = 0; i < THREADS; i++) {
            nn_ind[i] = split_hlen*i;
        }
        mpf_set_ui(tmp[0], 0);
        for (int i = 1; i < hlen; i++) {
            mpf_add(tmp[0], tmp[0], n[2*i]);
            mpf_mul(tmp[0], tmp[0], nn[i]);
        }
        mpf_add_ui(tmp[0], tmp[0], a+len+1);
        bfact_done++;
    }

    return NULL;
}

void *nextInd(int *index){
    *index = (*index + 1) % queueSize;
}

void *factorial(){
    long k = 2; 
    time_t raw_time = time(NULL);
    printf("\nStart time: %s",ctime(&raw_time));
    printf("----Start----\n");
    while(!done_calc){
        if(fullQueue && qfront == qrear) continue;
        fullQueue = 0;
        for(int i = 0; i < THREADS; i++){
            args[i] = (long*) malloc(sizeof(long) * 4);
            args[i][0] = k + split_len * i - 2;
            args[i][1] = k + split_len * (i + 1) - 2;
            args[i][2] = i;
            args[i][3] = k;
            pthread_create(&tid, NULL, start_bfact, (void *)(args[i]));
        }
        pthread_join(tid, NULL);
        while(bfact_done != THREADS){} //Wait for the last remaining thread to complete the factorial calculation
        bfact_done = 0;

        for(int i = 0; i < THREADS; i++) {
            threads_done[i] = 0;
        }
        pointerSwap(&tmpQueue[qfront], &tmp[0]);
        pointerSwap(&bfrQueue[qfront], &bfact_res[THREADS-1]);
        
        nextInd(&qfront);
        k += len;
        if(qfront == qrear) fullQueue = 1;
    }
}

void *approximate(){
    long iter = 0;
    int times = 1;
    while(!done_calc){
        if(!fullQueue && qrear == qfront) continue; //When the queue is empty
        mpf_div(tmpQueue[qrear], tmpQueue[qrear], bfrQueue[qrear]);
        mpf_add(eApprox, eApprox, tmpQueue[qrear]);
        if(mpf_cmp_d(tmpQueue[qrear], 0) == 0) done_calc = 1;

        nextInd(&qrear); 
        iter += len;
        if(iter > 100000*times){
            printf("iterations: %ld\n", iter);
            times++;
        }
        if(times >= 20) done_calc = 1;
    }
    t = clock() - t;
}

void *calculation(){
    pthread_t ctid;
    t = clock();
    pthread_create(&ctid, NULL, factorial, (void *) NULL);
    pthread_create(&ctid, NULL, approximate, (void *) NULL);
    while(!done_calc) sleep(2);
}

int main(int argc, char const *argv[]){ //arg 1: thread count
        // Core variables
    len = (int) pow(2.0, (double)LEN_POW);
    hlen = len >> 1;
    split_len = len / THREADS;
    split_hlen = split_len >> 1;
    ResultPrecison = (long) (log2(10) * pow(10, EXPO));
    FactorialPrecision = (long) (log2(10) * 300000);

        // Allocating memory
    args = (long **) malloc(sizeof(long*) * THREADS);
    n = (mpf_t*) malloc(mpf_t_size * len);
    nn = (mpf_t*) malloc(mpf_t_size * hlen); //n(x)*n(x+1) array
    nn_ind = (long *) malloc(sizeof(long) * THREADS); //An array of indexes of for threads
    bfact_res = (mpf_t*) malloc(mpf_t_size * THREADS); //Storing the results from the threads
    fact = (mpf_t*) malloc(mpf_t_size * THREADS);
    threads_done = (int *) malloc(sizeof(int) * THREADS);
    tmpQueue = (mpf_t*) malloc(mpf_t_size * queueSize);
    bfrQueue = (mpf_t*) malloc(mpf_t_size * queueSize);
    tmp = (mpf_t*) malloc(mpf_t_size);

        //Initializing variables
    mpf_init2(tmp[0], FactorialPrecision); //It's more than enough but i don't know how to calculate it's precision 
    for(int i = 0; i < queueSize; i++) {
        mpf_init2(tmpQueue[i], FactorialPrecision);
        mpf_init2(bfrQueue[i], FactorialPrecision);
    }
    for(int i = 0; i < len; i++) mpf_init(n[i]);
    for(int i = 0; i < hlen; i++) mpf_init2(nn[i], FactorialPrecision);
    for(int i = 0; i < THREADS; i++){
        nn_ind[i] = split_hlen*i; //index offset of length of the threads split
        threads_done[i] = 0;
        mpf_init2(bfact_res[i], FactorialPrecision);
        mpf_init2(fact[i], FactorialPrecision);
    }
    mpf_init2(eApprox, ResultPrecison);
    mpf_set_ui(eApprox, 1);
    printf("----Initial Variables----\n length power: %d\n base: %f\n exponent: %d\n approximate digits: %d\n threads: %d\n", LEN_POW, BASE, EXPO, PRECISION, THREADS);
    calculation();

        //Perfomance
    int mills = t * 1000 / CLOCKS_PER_SEC;
    int sec = mills/1000;
    int min = sec / 60;
    printf("----Finished----\n");
    printf("Time taken %d:%d:%d (mm:ss:ms)\n", min, sec%60, mills%1000);

        //Store result
    FILE *fp;
    fp = fopen("result.txt", "w");
    mpf_out_str(fp, 10, 0, eApprox);
    fseek(fp, 0, SEEK_SET);
    fputs("  ", fp);
    fclose(fp);

        //Free memory
    free(n);
    free(nn);
    free(nn_ind);
    free(bfact_res);
    free(fact);
    free(args);
    free(tmp);
    free(tmpQueue);
    free(bfrQueue);
    mpf_clear(eApprox);

        //Verify digits
    float million_digits = PRECISION / 1000000;
    if(million_digits >= 19.9998)
        system("\"C:/Program Files/nodejs/node\" ../Js/E/check.js ./result.txt ../e300M.txt");
    else if(million_digits >= 9.9998)
        system("\"C:/Program Files/nodejs/node\" ../Js/E/check.js ./result.txt ../e100M.txt");
    else
        system("\"C:/Program Files/nodejs/node\" ../Js/E/check.js ./result.txt ../e10M.txt");
    return 0;
}
