#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <unistd.h>
#include <pthread.h>

/* commands
gcc calculateE11.c -I/c/msys64/gmp -L/c/msys64/gmp/lib -lgmp -lpthread -o calc11
*/

#define EXPO 7
#define BASE 10.05
#define PRECISION (int)pow(BASE, EXPO)
#define THREADS 8
#define LEN_POW 16

int bsf_done = 0, final_thread_mul = 0, done_calc = 0, *threads_done, queueSize = 20,
    qfront = 0, qrear = 0, fullQueue = 0, queueByteSize = 0,
    qfront2 = 0, qrear2 = 0, fullQueue2 = 0, queueByteSize2 = 0,
    qfront3 = 0, qrear3 = 0, fullQueue3 = 0, queueByteSize3 = 0;

long len, hlen, *nn_ind, split_len, split_hlen, **args;
mpz_t *n, *nn, eApprox, prec, numerator,
    *bsfRes, *fact, *numQueue, *bsfQueue, *resQueue;
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
        mpz_init_set(temp_res, bsfRes[thread_ind]);
        bfact(m, b, thread_ind);
        mpz_mul(bsfRes[thread_ind], bsfRes[thread_ind], temp_res);

        mpz_clear(temp_res);
    }
    else
    {
        mpz_set_ui(fact[thread_ind], 1);
        for (long i = a + 1; i <= b; i++)
            mpz_mul_ui(fact[thread_ind], fact[thread_ind], i);
        mpz_set(bsfRes[thread_ind], fact[thread_ind]); // bsfRes is also being reseted with no overlap here
        mpz_set(nn[nn_ind[thread_ind]], fact[thread_ind]);
        mpz_set_ui(fact[thread_ind], 1);

        nn_ind[thread_ind]++;
    }
    return NULL;
}

void *start_bfact(void *vargv)
{ // Binary splitting factorial on thread level
    long *args = (long *)vargv;
    long a = args[0];   // Start index
    long b = args[1];   // End index
    long ind = args[2]; // Thread index
    long k = args[3];

    bfact(a, b, ind);

    int sibling = 1;
    if (ind != 0)
        bsf_done++;
    if ((ind + 1) % 2 == 0)
    { // Continues the BSF on the thread level like a tree from the bottom up
        // Each odd indexed thread multiplies it's values with his left sibling's value
        // recursively until there is only 1 thread remaining.
        // ie:(t0*t1),(t2*t3),(t4*t5),(t6*t7) -> (t1*t3), (t5*t7) -> (t3*t7) -> return t7
        while (1)
        {
            if (threads_done[ind - sibling])
            { // when the left sibling thread is done
                mpz_mul(bsfRes[ind], bsfRes[ind], bsfRes[ind - sibling]);
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
            mpz_set_ui(n[i], k + i); // pk - 1 + (0...p)
        while (bsf_done != THREADS - 1)
        {
        }
        for (int i = 0; i < THREADS; i++)
        {
            nn_ind[i] = split_hlen * i; // Sets next index of BSFactorial for each thread
        }
        for (int i = 1; i < hlen; i++)
        { // Calculates sequencial numerator
            // mpz_out_str(stdout, 10, n[2 * i]);
            // printf("\n");
            mpz_add(numerator, numerator, n[2 * i]);
            mpz_mul(numerator, numerator, nn[i]);
        }
        mpz_add_ui(numerator, numerator, a + len + 1);
        bsf_done++;
    }

    return NULL;
}

void *nextInd(int *index)
{
    *index = (*index + 1) % queueSize;
}

void *factorial()
{
    long k = 0;
    time_t raw_time = time(NULL);
    printf("\nStart time: %s", ctime(&raw_time));
    printf("----Start----\n");
    while (!done_calc)
    {
        if (fullQueue && qfront == qrear)
            continue;
        if (fullQueue2 && qfront2 == qrear2)
            continue;
        fullQueue = 0;
        fullQueue2 = 0;
        for (int i = 0; i < THREADS; i++)
        {
            args[i] = (long *)malloc(sizeof(long) * 4);
            args[i][0] = k + split_len * i;
            args[i][1] = k + split_len * (i + 1);
            args[i][2] = i;
            args[i][3] = k;
            pthread_create(&tid, NULL, start_bfact, (void *)(args[i]));
        }
        pthread_join(tid, NULL);
        while (bsf_done != THREADS)
        {
        } // Wait for the last remaining thread to complete the factorial calculation
        bsf_done = 0;

        for (int i = 0; i < THREADS; i++)
        {
            threads_done[i] = 0;
        }
        mpz_swap(numQueue[qfront2], numerator);
        mpz_swap(bsfQueue[qfront], bsfRes[THREADS - 1]);
        mpz_set_ui(numerator, 0);
        mpz_set_ui(bsfRes[THREADS - 1], 0);

        nextInd(&qfront2);
        nextInd(&qfront);
        k += len;
        if (qfront == qrear)
            fullQueue = 1;
        if (qfront2 == qrear2)
            fullQueue2 = 1;
    }
}

void *approximate()
{
    while (!done_calc)
    {
        if (!fullQueue && qrear == qfront)
            continue; // When the queue is empty
        mpz_fdiv_q(prec, prec, bsfQueue[qrear]);
        if (mpz_cmp_d(prec, 0) == 0)
        {
            done_calc = 1;
        }
        mpz_set(resQueue[qfront3], prec);
        fullQueue3 = 0;
        nextInd(&qfront3);

        mpz_set_ui(bsfQueue[qrear], 0);
        nextInd(&qrear);
    }
}

void *addApproximation()
{
    mpz_t res;
    long iter = 0;
    int times = 1;
    mpz_init(res);
    while (!done_calc)
    {
        if (!fullQueue2 && qrear2 == qfront2)
            continue; // When the queue is empty
        if (!fullQueue3 && qrear3 == qfront3)
            continue; // When the queue is empty
        mpz_mul(res, resQueue[qrear3], numQueue[qrear2]);
        mpz_add(eApprox, eApprox, res);

        mpz_set_ui(numQueue[qrear2], 0);
        mpz_set_ui(resQueue[qrear3], 0);
        nextInd(&qrear2);
        nextInd(&qrear3);

        iter += len;
        if (iter > 100000 * times)
        {
            printf("iterations: %ld\n", iter);
            times++;
        }
    }
    t = clock() - t;
    mpz_clear(res);
}

void *calculation()
{
    pthread_t ctid;
    t = clock();
    pthread_create(&ctid, NULL, factorial, (void *)NULL);
    pthread_create(&ctid, NULL, approximate, (void *)NULL);
    pthread_create(&ctid, NULL, addApproximation, (void *)NULL);
    while (!done_calc)
        sleep(2);
}

int main(int argc, char const *argv[])
{
    // Core variables
    len = (int)pow(2.0, (double)LEN_POW);
    hlen = len >> 1;
    split_len = len / THREADS;
    split_hlen = split_len >> 1;

    // Allocating memory
    args = (long **)malloc(sizeof(long *) * THREADS);
    n = (mpz_t *)malloc(mpz_t_size * len);
    nn = (mpz_t *)malloc(mpz_t_size * hlen);         // n(x)*n(x+1) array
    nn_ind = (long *)malloc(sizeof(long) * THREADS); // An array of indexes of for threads
    bsfRes = (mpz_t *)malloc(mpz_t_size * THREADS);  // Storing the results from the threads
    fact = (mpz_t *)malloc(mpz_t_size * THREADS);
    threads_done = (int *)malloc(sizeof(int) * THREADS);
    numQueue = (mpz_t *)malloc(mpz_t_size * queueSize);
    bsfQueue = (mpz_t *)malloc(mpz_t_size * queueSize);
    resQueue = (mpz_t *)malloc(mpz_t_size * queueSize);

    // Initializing variables
    mpz_init_set_ui(numerator, 0);
    for (int i = 0; i < queueSize; i++)
        mpz_inits(numQueue[i], bsfQueue[i], resQueue[i], NULL);
    for (int i = 0; i < len; i++)
        mpz_init(n[i]);
    for (int i = 0; i < hlen; i++)
        mpz_init(nn[i]);
    for (int i = 0; i < THREADS; i++)
    {
        nn_ind[i] = split_hlen * i; // index offset of length of the threads split
        threads_done[i] = 0;
        mpz_inits(bsfRes[i], fact[i], NULL);
    }
    mpz_init_set_d(prec, 10);
    mpz_pow_ui(prec, prec, PRECISION >> 1);
    mpz_pow_ui(prec, prec, 2);
    mpz_init_set(eApprox, prec);

    printf("----Initial Variables----\n length power: %d\n base: %f\n exponent: %d\n approximate digits: %d\n threads: %d\n", LEN_POW, BASE, EXPO, PRECISION, THREADS);
    calculation();

    // Perfomance
    int mills = t * 1000 / CLOCKS_PER_SEC;
    int sec = mills / 1000;
    int min = sec / 60;
    printf("----Finished----\n");
    printf("Time taken %d:%d:%d (mm:ss:ms)\n", min, sec % 60, mills % 1000);

    // Store result
    FILE *fp;
    fp = fopen("result.txt", "w");
    fwrite("a\n", 1, 2, fp);
    mpz_out_str(fp, 10, eApprox);
    fclose(fp);
    printf("Stored the constant in result.txt\n");

    // Free memory
    for (int i = 0; i < len; i++)
        mpz_clear(n[i]);
    for (int i = 0; i < hlen; i++)
        mpz_clear(nn[i]);
    free(n);
    free(nn);
    free(nn_ind);
    for (int i = 0; i < THREADS; i++)
        mpz_clears(bsfRes[i], fact[i], NULL);
    free(bsfRes);
    free(fact);
    free(args);
    for (int i = 0; i < queueSize; i++)
        mpz_clears(numQueue[i], bsfQueue[i], resQueue[i], NULL);
    free(numQueue);
    free(bsfQueue);
    free(resQueue);
    mpz_clears(prec, eApprox, numerator, NULL);

    // Verify digits
    float million_digits = PRECISION / 1000000;
    if (million_digits <= 10)
        system("\"C:/Program Files/nodejs/node\" ../Js/E/check.js ./result.txt ../e10M.txt");
    else if (million_digits <= 100)
        system("\"C:/Program Files/nodejs/node\" ../Js/E/check.js ./result.txt ../e100M.txt");
    else if (million_digits <= 300)
        system("\"C:/Program Files/nodejs/node\" ../Js/E/check.js ./result.txt ../e300M.txt");
    else
        system("\"C:/Program Files/nodejs/node\" ../Js/E/check.js ./result.txt ../e500M.txt");
    return 0;
}
