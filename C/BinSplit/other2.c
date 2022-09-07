#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <unistd.h>

mpz_t prod, q_f_res, sum, p_f_res;

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

void* P_F(int a, int b) {
	if (!((a + b) % 2) && b - a > 2) {
        int m = (a + b) / 2;
		mpz_t temp_res;

		P_F(a, m);
		mpz_init_set(temp_res, p_f_res);
		mpz_set_ui(q_f_res, 1);
		Q_F(m, b);
		mpz_mul(temp_res, temp_res, q_f_res);
		P_F(m, b);
		mpz_add(p_f_res, temp_res, p_f_res);

		mpz_clear(temp_res);
	} else {
		mpz_set_ui(sum, 0);
		for (int i = a + 1; i <= b; i++){
			mpz_set_ui(q_f_res, 1);
			Q_F(i, b);
			mpz_add(sum, sum, q_f_res);
		}
		mpz_set(p_f_res, sum);
	}
}

int main(int argc, char const *argv[]){
	long int digits = 100000;
    const int k = digits / 4;
	printf("k:%d\n", k);
    int a = 0;
    int b = k;

    mpz_t Q, P, E, exponent, base;
	mpz_init(Q);
	mpz_init(P);
	mpz_init(E);
	mpz_init(prod);
	mpz_init(sum);
	mpz_init(exponent);
	mpz_init_set_ui(base, 10);
	mpz_pow_ui(exponent, base, digits);
	mpz_init_set_ui(q_f_res, 1);
	mpz_init_set_ui(p_f_res, 1); //maybe set 0? 

		//Calculate Q(a,b)
	clock_t start = clock();
	Q_F(a, b);
	mpz_set(Q, q_f_res);
		//Calculate P(a,b)
	P_F(a, b);
	mpz_mul(P, p_f_res, exponent);
		//Calculate E
	mpz_div(E, P, Q);
	printf("Run Time: %dms\n", clock() - start);
		//Output to file
	FILE *fp;
	fp = fopen("result.txt", "w+");
	if(a) fprintf(fp, "2");
	mpz_out_str(fp, 10, E);

	fclose(fp);

	mpz_clears(Q, P, E, prod, sum,q_f_res, p_f_res, exponent, base, NULL);
	system("\"C:/Program Files/nodejs/node\" ../../Js/E/check.js ./result.txt ../../e10M.txt");

    return 0;
}
