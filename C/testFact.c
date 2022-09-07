#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <unistd.h>
#include <pthread.h>

/*
gcc testFact.c -I/c/msys64/gmp -L/c/msys64/gmp/lib -lgmp -lpthread -o test
*/
#define EXPO 8
#define BASE 11.52
#define PRECISION (int)pow(BASE, EXPO)
#define THREADS 8
#define LEN_POW 21

int bfact_done = 0, final_thread_mul = 0, done_calc = 0, *threads_done,
    queueSize = 20, qfront = 0, qrear = 0, fullQueue = 0;
long len, hlen, *nn_ind, split_len, split_hlen, **args;
mpz_t *n, *nn, eApprox, prec, *tmp,
    *bfact_res, *fact, *tmpQueue, *bfrQueue;
pthread_t tid;
clock_t t;

int mpz_t_size = sizeof(mpz_t);

void *bfact(long a, long b, long thread_ind)
{ // Binary splitting factorial
    if (b - a > 2)
    {
        long m = (long)floor((a + b) / 2);
        mpz_t temp_res;

        bfact(a, m, thread_ind);
        mpz_init_set(temp_res, bfact_res[thread_ind]);
        bfact(m, b, thread_ind);
        mpz_mul(bfact_res[thread_ind], bfact_res[thread_ind], temp_res);

        mpz_clear(temp_res);
    }
    else
    {
        mpz_set_ui(fact[thread_ind], 1);
        for (long i = a + 1; i <= b; i++)
            mpz_mul_ui(fact[thread_ind], fact[thread_ind], i);
        mpz_set(bfact_res[thread_ind], fact[thread_ind]); // bfact_res is also being reseted with no overlap here
        mpz_set(nn[nn_ind[thread_ind]], fact[thread_ind]);
        nn_ind[thread_ind]++;
    }
    return NULL;
}

void *start_bfact(void *vargv)
{
    long *args = (long *)vargv;
    long a = args[0];
    long b = args[1];
    long ind = args[2];
    long k = args[3];

    bfact(a, b, ind);

    int sibling = 1;
    if (ind != 0)
        bfact_done++;
    if ((ind + 1) % 2 == 0)
    {
        // Each odd indexed thread multiplies it's values with his left sibling's value
        // recursively until there is only 1 thread remaining.
        // ie:(t0*t1),(t2*t3),(t4*t5),(t6*t7) -> (t1*t3), (t5*t7) -> (t3*t7) -> return t7
        while (1)
        {
            if (threads_done[ind - sibling])
            { // if left sibling thread is done
                mpz_mul(bfact_res[ind], bfact_res[ind], bfact_res[ind - sibling]);
                sibling *= 2;
                threads_done[ind] = 1; // prevent overlapping
                if (((ind + 1) / sibling) % 2 == 1)
                    break;
                threads_done[ind] = 0;
            }
        }
    }
    else
        threads_done[ind] = 1;

    if (ind == 0)
    {
        // Make the first thread
        for (int i = 0; i < len; i++)
            mpz_set_d(n[i], k - 1 + i);
        while (bfact_done != THREADS - 1)
        {
        }
        for (int i = 0; i < THREADS; i++)
        {
            nn_ind[i] = split_hlen * i;
        }
        mpz_set_d(tmp[0], 0);
        for (int i = 1; i < hlen; i++)
        {
            mpz_add(tmp[0], tmp[0], n[2 * i]);
            mpz_mul(tmp[0], tmp[0], nn[i]);
        }
        mpz_add_ui(tmp[0], tmp[0], a + len + 1);
        bfact_done++;
    }

    return NULL;
}

void *factorial()
{
    long k = 2;
    time_t raw_time = time(NULL);
    printf("\nStart time: %s", ctime(&raw_time));
    printf("----Start----\n");
    t = clock();
    while (!done_calc)
    {
        for (int i = 0; i < THREADS; i++)
        {
            args[i] = (long *)malloc(sizeof(long) * 4);
            args[i][0] = k + split_len * i - 2;
            args[i][1] = k + split_len * (i + 1) - 2;
            args[i][2] = i;
            args[i][3] = k;
            pthread_create(&tid, NULL, start_bfact, (void *)(args[i]));
        }
        pthread_join(tid, NULL);
        while (bfact_done != THREADS)
        {
        } // Wait for the last remaining thread to complete the factorial calculation
        bfact_done = 0;
        done_calc = 1;
    }
}

void *calculation()
{
    pthread_t ctid;
    pthread_create(&ctid, NULL, factorial, (void *)NULL);
    while (!done_calc)
        sleep(1);
}

int main(int argc, char const *argv[])
{ // arg 1: thread count
    // Core variables
    len = (int)pow(2.0, (double)LEN_POW);
    hlen = len >> 1;
    split_len = len / THREADS;
    split_hlen = split_len >> 1;

    // Allocating memory
    args = (long **)malloc(sizeof(long *) * THREADS);
    n = (mpz_t *)malloc(mpz_t_size * len);
    nn = (mpz_t *)malloc(mpz_t_size * hlen);           // n(x)*n(x+1) array
    nn_ind = (long *)malloc(sizeof(long) * THREADS);   // An array of indexes of for threads
    bfact_res = (mpz_t *)malloc(mpz_t_size * THREADS); // Storing the results from the threads
    fact = (mpz_t *)malloc(mpz_t_size * THREADS);
    threads_done = (int *)malloc(sizeof(int) * THREADS);
    tmpQueue = (mpz_t *)malloc(mpz_t_size * queueSize);
    bfrQueue = (mpz_t *)malloc(mpz_t_size * queueSize);
    tmp = (mpz_t *)malloc(mpz_t_size);

    // Initializing variables
    mpz_init(tmp[0]);
    for (int i = 0; i < queueSize; i++)
        mpz_inits(tmpQueue[i], bfrQueue[i], NULL);
    for (int i = 0; i < len; i++)
        mpz_init(n[i]);
    for (int i = 0; i < hlen; i++)
        mpz_init(nn[i]);
    for (int i = 0; i < THREADS; i++)
    {
        nn_ind[i] = split_hlen * i; // index offset of length of the threads split
        threads_done[i] = 0;
        mpz_inits(bfact_res[i], fact[i], NULL);
    }

    printf("----Initial Variables----\n length power: %d\n base: %f\n exponent: %d\n approximate digits: %d\n threads: %d\n", LEN_POW, BASE, EXPO, PRECISION, THREADS);
    calculation();

    // Perfomance
    int mills = t * 1000 / CLOCKS_PER_SEC;
    int sec = mills / 1000;
    int min = sec / 60;
    printf("----Finished----\n");
    printf("Time taken %d:%d:%d (mm:ss:ms)\n", min, sec % 60, mills % 1000);

    // Free memory
    free(n);
    free(nn);
    free(nn_ind);
    free(bfact_res);
    free(fact);
    free(args);
    free(tmp);
    free(tmpQueue);
    free(bfrQueue);

    return 0;
}
