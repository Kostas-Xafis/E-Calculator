//!BigInt serialization
const toJson = data => JSON.stringify(data, (key, value) => (typeof value === "bigint" ? value.toString() + "n" : value), 4);

const toObject = json => {
	JSON.parse(json, (key, value) => {
		if (typeof value === "string" && /^\d+n$/.test(value)) {
			return BigInt(value.substr(0, value.length - 1));
		}
		return value;
	});
};

module.exports = { toJson, toObject };
