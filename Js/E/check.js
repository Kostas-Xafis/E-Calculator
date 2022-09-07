const fs = require("fs/promises");
const { open } = fs;
const { argv } = require("process");

(async function () {
	let bytes = 1024 * 1024 * 10; // 10MB
	// const fd1 = await open("../../c/result.txt", "r");
	// const fd2 = await open("../../e300M.txt", "r");
	const fd1 = await open(argv[2], "r");
	const fd2 = await open(argv[3], "r");
	let i = 0;
	while (true) {
		const buffer = Buffer.alloc(bytes);
		const buf1 = await fd1.read({ buffer, position: bytes * i });
		const buf2 = await fd2.read({ buffer, position: bytes * i });
		const str1 = buf1.buffer.toString();
		const str2 = buf2.buffer.toString();
		if (buf1.bytesRead != bytes && i == 0) {
			bytes >>= 10;
			continue;
		}
		if (str1 != str2 || buf1.bytesRead != bytes) break;
		i++;
	}
	console.log(`Total correct digits: ${i * bytes}`);
})();
