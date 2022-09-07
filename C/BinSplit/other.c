#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

long long Q_F(int a, int b) {
	// printf("a:%d\tb:%d\n", a, b);
	if (!((a + b) % 2) && b - a > 2) {
        int m = (int) (a + b) / 2;
		return  Q_F(a, m) * Q_F(m, b);
	} else {
		long long prod = 1;
		for (int i = a + 1; i <= b; i++) prod = prod * i;
		// printf("prod:%d\n", prod);
		return prod;
	}
}

long long P_F(int a, int b) {
	// printf("a:%d\tb:%d\n", a, b);
	if (!((a + b) % 2) && b - a > 2) {
        int m = (a + b) / 2;
		return P_F(a, m) * Q_F(m, b) + P_F(m, b);
	} else {
		long long sum = 0;
		for (int i = a + 1; i <= b; i++) sum += Q_F(i, b);
		// printf("sum:%ld\n", sum);
		return sum;
	}
}

int main(int argc, char const *argv[]){
    const int n = atoi(argv[argc-1]);
    int a = n % 2;
    int b = n;

    long long Q = Q_F(a, b);
    long long P = P_F(a, b);
    printf("Q:%lld\t P:%lld\n", Q, P);

    double E = ((double) P / (double) Q) + 1 + a;
    printf("e:%.20lf\n", E);
    
	

    return 0;
}
