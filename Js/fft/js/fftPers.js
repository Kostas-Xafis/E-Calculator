//?Bitwise OR (|) seems to be on of few ways to cast a float
//? to an int, that also covers edge cases like 1e-10 or
const { appendFile } = require("fs");
const { performance } = require("perf_hooks");
const Complex = require("./complex.js");
const { add, sub, mul, pow } = Complex;

//! The accuracy must be fine tuned
const ROUA = 1000_000_000; //Roots Of Unity Accuracy
const b_ROUA = BigInt(ROUA); //Roots Of Unity Accuracy
Complex.setAccuracy(b_ROUA);

const RootsOfUnity = (fft_len => {
	function unity_roots(n, inverse = 1) {
		const t = (inverse * 2 * Math.PI) / n;
		const [cos_r, sin_r] = [Math.cos(t), Math.sin(t)];
		const base = new Complex({
			re: cos_r | 0,
			im: sin_r | 0,
			re_rem: (cos_r * ROUA) % ROUA | 0,
			im_rem: (sin_r * ROUA) % ROUA | 0
		});
		return Array(n)
			.fill(0) // Assigning omega to array
			.map((x, i) => pow(base.clone(), i)); // raise omega to i-th power
	}

	let len = fft_len,
		roots = {},
		inv_roots = {};
	while (len >= 1) {
		roots[len] = unity_roots(len);
		inv_roots[-len] = unity_roots(len, -1); //Roots for the IFFT
		len >>= 1;
	}

	return { ...roots, ...inv_roots };
})(2048);
console.log(RootsOfUnity);
let inverse = 1;
let trues = 0,
	outof = 0;
function fft(P) {
	const n = P.length;
	if (n == 1) return P;

	const omArr = RootsOfUnity[n * inverse]; // omega array
	const [Pe, Po] = [P.filter((_x, i) => !(i % 2)), P.filter((_x, i) => i % 2)]; //Split evens-odds
	const [ye, yo] = [fft(Pe), fft(Po)];
	let y = Array(n).fill(0);

	for (let i = 0; i < n / 2; i++) {
		y[i] = add(ye[i], mul(omArr[i], yo[i]));
		y[i + n / 2] = sub(ye[i], mul(omArr[i], yo[i]));
	}
	return y;
}

//!Start Here
// for (let k = 1; k <= 49; k++) {
// 	for (let j = k + 1; j <= 50; j++) {
// 		inverse = 1;
// 		let [a, b] = [BigInt(k) ** 250n, BigInt(j) ** 250n];
// 		const real_result = a * b;
// 		let m1 = [];
// 		let m2 = [];

// 		while (a >= 1n || b >= 1n) {
// 			m1.push(new Complex({ re: a % 10n, im: 0n }));
// 			a /= 10n;
// 			m2.push(new Complex({ re: b % 10n, im: 0n }));
// 			b /= 10n;
// 		}

// 		const [fft_len, fft_len_pow2] = (() => {
// 			let i = 1;
// 			let pow2 = 1;
// 			while ((i <<= 1) < 2 * m1.length) pow2++;
// 			return [i, pow2];
// 		})(); //fft_size must be a power of 2 and bigger than at least 2*max(d) (d:digits of a\b)
// 		const padding = fft_len - m1.length;
// 		m1 = m1.concat(Array(padding).fill(new Complex({ re: 0n, im: 0n })));
// 		m2 = m2.concat(Array(padding).fill(new Complex({ re: 0n, im: 0n })));

// 		let start = performance.now();
// 		const f1 = fft(m1);
// 		const f2 = fft(m2);
// 		const f_mul = f1.map((x, i) => mul(x, f2[i]));
// 		inverse = -1; // Changing FFT algorithm to the IFFT
// 		const fft_res = fft(f_mul);

// 		//! Removing rounding errors (very important)
// 		const per90 = (Complex.acc * 60n) / 100n; //!Also the percentage must be fine tuned
// 		for (let i = 0; i < fft_res.length; i++) {
// 			if (fft_res[i].re_rem >= per90) fft_res[i].re += 1n;
// 		}

// 		let sum = 0n;
// 		for (let i = 0n; i < fft_res.length; i++) sum += fft_res[i].re * 10n ** i;

// 		const result = sum >> BigInt(fft_len_pow2);
// 		let correct = result === real_result;
// 		trues += correct === true;
// 		console.log({ time: Math.floor(performance.now() - start), j, k, correct });
// 	}
// 	outof += k;
// }
appendFile("correct.txt", `\nCorrect:${trues}\t out of:${outof}`, console.error);
