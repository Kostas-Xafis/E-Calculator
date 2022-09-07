// void *m8(long n, mpz_t* res){
//     mpz_set_ui(nm1, n - 1);
//     mpz_set_ui(n7, n+ 7);

//     //  n1n2 = (n + 1n) * (n + 2n);
//     mpz_set_ui(n1, n + 1);
//     mpz_set_ui(n2, n+ 2);
// 	mpz_mul(n1n2, n2, n1);

//     //n3n4 = (n + 3n) * (n + 4n);
//     mpz_set_ui(n3, n+ 3);
//     mpz_set_ui(n4, n+ 4);
// 	mpz_mul(n3n4, n3, n4);

// 	//n5n6 = (n + 5n) * (n + 6n);
//     mpz_set_ui(n5, n +5);
//     mpz_set_ui(n6, n+ 6);
// 	mpz_mul(n5n6, n5, n6);
	
//     //n3n4n5n6 = n3n4 * n5n6;
// 	mpz_mul(n3n4n5n6, n3n4, n5n6);
    
// 	// acc /= (n - 1n) * n * n1n2 * n3n4n5n6;
// 	mpz_mul_ui(tmp1, nm1, n);

// 	mpz_mul(tmp2, n1n2, n3n4n5n6);
//     mpz_mul(tmp2, tmp2, tmp1);
//     mpz_div(prec, prec, tmp2);

//     // return acc * ((n1n2 * (n + 1n) + n + 3n) * n3n4n5n6 + (n5n6 * (n + 5n) + n + 7n));

//     mpz_mul(tmp1, n1n2, n1);
//     mpz_add(tmp1, tmp1, n3);
//     mpz_mul(tmp1, tmp1, n3n4n5n6);

//     mpz_mul(tmp2, n5n6, n5);
//     mpz_add(tmp2, tmp2, n7);

//     mpz_add(tmp2, tmp2, tmp1);
//     mpz_mul(*res, tmp2, prec);
//     return NULL;
// }