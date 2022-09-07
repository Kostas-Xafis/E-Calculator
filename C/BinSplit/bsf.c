#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <unistd.h>
//	gcc bsf.c -lgmp -o bsf

//Demonstrating the binary splitting algorithm in factorial function
mpz_t prod, q_f_res;

void* Q_F(int a, int b) {
	if (!((a + b) % 2) && b - a > 2) {
        int m = (int) (a + b) / 2;
		mpz_t temp_res;
		Q_F(a, m);
		mpz_init_set(temp_res, q_f_res);
		Q_F(m, b);
		mpz_mul(q_f_res, q_f_res, temp_res);
		mpz_clear(temp_res);
	} else {
		mpz_set_ui(prod, 1);
		for (int i = a + 1; i <= b; i++) mpz_mul_ui(prod, prod, i);
		mpz_set(q_f_res, prod);
	}
}

void* fact(int n){
	while(n-- > 2)
		mpz_mul_ui(prod, prod, n+1);
}

int main(int argc, char const *argv[]){
    const int n = atoi(argv[argc-1]);;
    int a = n % 2;
    int b = n;

    mpz_t Q;
	mpz_init(Q);
	mpz_init(prod);
	mpz_init_set_ui(q_f_res, 1);

	clock_t start = clock();
	Q_F(a, b);
	printf("Binary fact took:%d ms\n", clock() - start);
	
    mpz_set(Q, q_f_res);
	// mpz_set_ui(prod, 1);

	start = clock();
	fact(n);
	printf("Normal fact took:%d ms\n", clock() - start);

	FILE *fp;
	fp = fopen("fact.txt", "w+");
	mpz_out_str(fp, 10, Q);

	fclose(fp);

	mpz_clears(Q, prod, q_f_res);    
    return 0;
}
