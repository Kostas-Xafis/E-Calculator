const n = 10000n,
	a = n % 2n,
	b = n;
function Q_F(a, b) {
	if (!((a + b) % 2n) && b - a > 2n) {
		return Q_F(a, (a + b) / 2n) * Q_F((a + b) / 2n, b);
	} else {
		let prod = 1n;
		for (let i = a + 1n; i <= b; i++) prod = prod * i;
		return prod;
	}
}

function P_F(a, b) {
	if (!((a + b) % 2n) && b - a > 2n) {
		let m = (a + b) / 2n;
		return P_F(a, m) * Q_F(m, b) + P_F(m, b);
	} else {
		let sum = 0n;
		for (let i = a + 1n; i <= b; i++) {
			sum += Q_F(i, b);
		}
		return sum;
	}
}

function fact(n) {
	if (n == 2n) return 2n;
	return n * fact(n - 1n);
}

// console.time("Regular Factorial");
// for (let i = 0; i < 10_000; i++) fact(n);
// console.timeEnd("Regular Factorial");

// let res = 0n;
// console.time("Binary Factorial");
// for (let i = 0; i < 10_000; i++) res = Q_F(a, b);
// console.timeEnd("Binary Factorial");

// console.log(res);

console.time("Binary Split 100x");
for (let i = 0; i < 10; i++) {
	const Q = Q_F(a, b);
	const P = P_F(a, b);
	const res = P / Q + 1n + a;
	console.log(i);
}
console.timeEnd("Binary Split 100x");
