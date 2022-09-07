const fs = require("fs");
const Complex = require("./complex");
const { toJson } = require("./bigintToJson");
const _sin = Complex.sin;
const _cos = Complex.cos;

let n = 2n ** 10n;
let result = {};
Complex.setAccuracy(1000000n);
const t1 = 2n * BigInt(Math.floor(Math.PI * 100000000));
while (n >= 1n) {
	let t = t1 / n;
	let [sin, cos] = [_sin(t), _cos(t)];
	result[n] = { sin, cos };
	result[-n] = { sin, cos };
	n >>= 1n;
}

fs.writeFile("SinCos.txt", toJson(result), { encoding: "utf-8" }, console.error);
