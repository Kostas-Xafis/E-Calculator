var _c;
module.exports = (_c = class Complex {
        constructor(n) {
            this.re = typeof n.re == "bigint" ? n.re : BigInt(n.re);
            if (Number.isNaN(n.re_rem) || n.re_rem === undefined)
                this.re_rem = 0n;
            else
                this.re_rem = typeof n.re_rem == "bigint" ? n.re_rem : BigInt(n.re_rem);
            this.im = typeof n.im == "bigint" ? n.im : BigInt(n.im);
            if (Number.isNaN(n.im_rem) || n.im_rem === undefined)
                this.im_rem = 0n;
            else
                this.im_rem = typeof n.im_rem == "bigint" ? n.im_rem : BigInt(n.im_rem);
        }
        static add(a, b) {
            const [_a, _b, c, d] = [
                a.re * Complex.acc + a.re_rem,
                a.im * Complex.acc + a.im_rem,
                b.re * Complex.acc + b.re_rem,
                b.im * Complex.acc + b.im_rem
            ];
            let re = _a + c;
            let im = _b + d;
            let re_rem = re % Complex.acc;
            let im_rem = im % Complex.acc;
            re /= Complex.acc;
            im /= Complex.acc;
            return new Complex({ re, im, re_rem, im_rem });
        }
        static sub(a, b) {
            const [_a, _b, c, d] = [
                a.re * Complex.acc + a.re_rem,
                a.im * Complex.acc + a.im_rem,
                b.re * Complex.acc + b.re_rem,
                b.im * Complex.acc + b.im_rem
            ];
            let re = _a - c;
            let im = _b - d;
            let re_rem = re % Complex.acc;
            let im_rem = im % Complex.acc;
            re /= Complex.acc;
            im /= Complex.acc;
            return new Complex({ re, im, re_rem, im_rem });
        }
        static mul(a, b) {
            const [_a, _b, c, d] = [
                a.re * Complex.acc + a.re_rem,
                a.im * Complex.acc + a.im_rem,
                b.re * Complex.acc + b.re_rem,
                b.im * Complex.acc + b.im_rem
            ];
            const ac_bd = _a * c - _b * d, ad_bc = _a * d + _b * c;
            const re = ac_bd / Complex.acc_pow2, im = ad_bc / Complex.acc_pow2, re_rem = (ac_bd % Complex.acc_pow2) / Complex.acc, im_rem = (ad_bc % Complex.acc_pow2) / Complex.acc;
            return new Complex({ re, im, re_rem, im_rem });
        }
        static sin(n) {
            let sum = 1n;
            let fact = 6n;
            let expo = n ** 3n;
            let sign = -1n;
            const np2 = n ** 2n;
            for (let i = 3n; i < 50n; i += 2n) {
                sum += sign * (expo / fact);
                expo *= np2;
                fact *= (i + 2n) * (i + 1n);
                sign *= -1n;
            }
            return sum;
        }
        static cos(n) {
            let sum = 1n;
            let fact = 2n;
            let expo = n ** 2n;
            let sign = -1n;
            const np2 = n ** 2n;
            for (let i = 2n; i < 50n; i += 2n) {
                sum += sign * (expo / fact);
                expo *= np2;
                fact *= (i + 2n) * (i + 1n);
                sign *= -1n;
            }
            console.log(sum);
            return sum;
        }
        static pow(a, b) {
            if (b == 0)
                return new Complex({ re: 1n, im: 0n, re_rem: 0n, im_rem: 0n });
            if (b == 1)
                return a.clone();
            let pow = a.clone();
            for (let i = 1; i < b; i++)
                pow = Complex.mul(pow, a);
            return pow;
        }
        clone() {
            return new Complex({
                re: this.re,
                im: this.im,
                re_rem: this.re_rem,
                im_rem: this.im_rem
            });
        }
        static setAccuracy(a) {
            Complex.acc = a;
            Complex.acc_pow2 = a ** 2n;
        }
    },
    _c.acc = 1n,
    _c.acc_pow2 = 1n,
    _c);
