const fs = require("fs"),
	cluster = require("cluster"),
	{ exec } = require("child_process"),
	exponent = 100_000,
	forks = 3,
	forks_b = BigInt(forks);
let workersHandler,
	messageHandler,
	done = 0,
	precision,
	eApproximation = 0n;

function findDigits(n) {
	if (n < 0) return 0;

	if (n <= 1) return 1;

	let x = n * Math.log10(n / Math.E) + Math.log10(2 * Math.PI * n) / 2.0;

	return Math.floor(x) + 1;
}

function main(iter, offset) {
	let eApprox = 0n,
		i = offset + 1n,
		prec = precision,
		times = 1n,
		subFact = 1n,
		k = 1n;
	for (; k <= i; subFact *= k++);

	for (; i <= iter; ) {
		prec /= subFact;
		eApprox += prec;

		subFact = 1n;
		for (k = i + 1n; k <= i + forks_b; k += 1n) subFact *= k;
		i += forks_b;

		if (prec == 0n) break;
		if (i >= 5000n * times) {
			console.log(5000n * times);
			times += 1n;
		}
	}
	return eApprox;
}

function results(res) {
	fs.writeFile("res.txt", String(res), { encoding: "utf-8" }, err => (err ? console.log(err) : console.log("Done")));
	exec("node ./check.js ./res.txt", (err, stdo, stderr) => {
		if (err) console.error(err);
		if (stdo) console.error(stdo);
		if (stderr) console.error(stderr);
	});
}

if (cluster.isMaster) {
	console.log(`Master ${process.pid} is running`);

	workersHandler = () => {
		messageHandler = res => {
			done++;
			eApproximation += BigInt(res.res);
			cluster.workers[res.id].kill();
			if (done == forks) {
				console.timeEnd("E");
				results(String(eApproximation));
			}
		};
		//?Fork workers
		for (let i = 0; i < forks; i++) {
			const worker = cluster.fork();
			worker.on("message", messageHandler);
			worker.send(i);
		}
		console.time("E");
	};
	workersHandler();
} else {
	//* Worker code
	precision = 10n ** BigInt(exponent);
	process.on("message", offset => {
		let result = main(BigInt(findDigits(exponent)), BigInt(offset));
		process.send({ res: String(result), id: cluster.worker.id });
	});
}
