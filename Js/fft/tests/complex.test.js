const Complex = require("../js/complex.js");
const { add, sub, mul, pow, div } = Complex;
Complex.setAccuracy(1000n);

//! ADDITION

test("Adds 2 BigInt Floating Complex Numbers #1", () => {
	let a = new Complex({ re: 5, re_rem: 500, im: -2, im_rem: -800 });
	let b = new Complex({ re: 1, re_rem: 700, im: 4, im_rem: 600 });
	let res = add(a, b);
	expect(res.re).toBe(7n);
	expect(res.re_rem).toBe(200n);
	expect(res.im).toBe(1n);
	expect(res.im_rem).toBe(800n);
});

test("Adds 2 BigInt Floating Complex Numbers #2", () => {
	let a = new Complex({ re: -5, re_rem: -500, im: 2, im_rem: 800 });
	let b = new Complex({ re: -2, re_rem: -700, im: -4, im_rem: -600 });
	let res = add(a, b);
	expect(res.re).toBe(-8n);
	expect(res.re_rem).toBe(-200n);
	expect(res.im).toBe(-1n);
	expect(res.im_rem).toBe(-800n);
});

//! SUBTRACTION

test("Subtracts 2 BigInt Floating Complex Numbers #1", () => {
	let a = new Complex({ re: 5, re_rem: 500, im: -2, im_rem: -800 });
	let b = new Complex({ re: 1, re_rem: 700, im: 4, im_rem: 600 });
	let res = sub(a, b);
	expect(res.re).toBe(3n);
	expect(res.re_rem).toBe(800n);
	expect(res.im).toBe(-7n);
	expect(res.im_rem).toBe(-400n);
});

test("Subtracts 2 BigInt Floating Complex Numbers #2", () => {
	let a = new Complex({ re: -5, re_rem: -500, im: 2, im_rem: 800 });
	let b = new Complex({ re: -2, re_rem: -700, im: -4, im_rem: -600 });
	let res = sub(a, b);
	expect(res.re).toBe(-2n);
	expect(res.re_rem).toBe(-800n);
	expect(res.im).toBe(7n);
	expect(res.im_rem).toBe(400n);
});

//! MULTIPLICATION

test("Multiplies 2 BigInt Floating Complex Numbers #1", () => {
	let a = new Complex({ re: 5, re_rem: 500, im: -2, im_rem: -800 });
	let b = new Complex({ re: 1, re_rem: 700, im: 4, im_rem: 600 });
	let res = mul(a, b);
	expect(res.re).toBe(22n);
	expect(res.re_rem).toBe(230n);
	expect(res.im).toBe(20n);
	expect(res.im_rem).toBe(540n);
});

test("Multiplies 2 BigInt Floating Complex Numbers #2", () => {
	let a = new Complex({ re: -5, re_rem: -500, im: 2, im_rem: 800 });
	let b = new Complex({ re: -2, re_rem: -700, im: -4, im_rem: -600 });
	let res = mul(a, b);
	expect(res.re).toBe(27n);
	expect(res.re_rem).toBe(730n);
	expect(res.im).toBe(17n);
	expect(res.im_rem).toBe(740n);
});
