using MultiPrecision;
using MultiPrecisionComplex;

namespace ComplexBessel {
    public class HankelExpansion<N> where N : struct, IConstant {
        public MultiPrecision<N> Nu { get; }

        private readonly List<MultiPrecision<N>> a_coef;

        private readonly MultiPrecision<N> cospi_nu, sinpi_nu;

        public HankelExpansion(MultiPrecision<N> nu) {
            Nu = nu;
            a_coef = [1];
            cospi_nu = MultiPrecision<N>.CosPI(Nu);
            sinpi_nu = MultiPrecision<N>.SinPI(Nu);
        }

        public MultiPrecision<N> ACoef(int n) {
            for (int k = a_coef.Count; k <= n; k++) {
                MultiPrecision<N> a = a_coef.Last() * (4 * Nu * Nu - checked((2 * k - 1) * (2 * k - 1))) / (k * 8);
                a_coef.Add(a);
            }

            return a_coef[n];
        }

        public Complex<N> Omega(Complex<N> z) {
            Complex<N> omega = z - MultiPrecision<N>.Ldexp(2 * Nu + 1, -2) * MultiPrecision<N>.PI;

            return omega;
        }

        public (Complex<N> c_even, Complex<N> c_odd) BesselJYCoef(Complex<N> z, int max_term = 256) {
            Complex<N> v = 1 / (z * z), w = -v;

            Complex<N> c_even = ACoef(0), c_odd = ACoef(1);

            for (int k = 1; k <= max_term; k++) {
                Complex<N> dc_even = w * ACoef(2 * k);
                Complex<N> dc_odd = w * ACoef(2 * k + 1);

                c_even += dc_even;
                c_odd += dc_odd;

                if ((c_even.Exponent - dc_even.Exponent >= MultiPrecision<N>.Bits) &&
                    (c_odd.Exponent - dc_odd.Exponent >= MultiPrecision<N>.Bits)) {

                    return (c_even, c_odd / z);
                }

                w *= -v;
            }

            return (MultiPrecision<N>.NaN, MultiPrecision<N>.NaN);
        }

        public Complex<N> BesselICoef(Complex<N> z, int max_term = 256) {
            Complex<N> v = 1 / z, w = -v;

            Complex<N> c = ACoef(0);

            for (int k = 1; k <= max_term; k++) {
                Complex<N> dc = w * ACoef(k);

                c += dc;

                if (c.Exponent - dc.Exponent >= MultiPrecision<N>.Bits) {
                    return c;
                }

                w *= -v;
            }

            return MultiPrecision<N>.NaN;
        }

        public Complex<N> BesselKCoef(Complex<N> z, int max_term = 256) {
            Complex<N> v = 1 / z, w = v;

            Complex<N> c = ACoef(0);

            for (int k = 1; k <= max_term; k++) {
                Complex<N> dc = w * ACoef(k);

                c += dc;

                if (c.Exponent - dc.Exponent >= MultiPrecision<N>.Bits) {
                    return c;
                }

                w *= v;
            }

            return MultiPrecision<N>.NaN;
        }

        public Complex<N> BesselJ(Complex<N> z) {
            if (z.R.Sign == Sign.Minus) {
                return (cospi_nu, (z.I.Sign == Sign.Plus) ? sinpi_nu : -sinpi_nu) * BesselJ(-z);
            }

            (Complex<N> c_even, Complex<N> c_odd) = BesselJYCoef(z);
            if (Complex<N>.IsNaN(c_even) || Complex<N>.IsNaN(c_odd)) {
                return Complex<N>.NaN;
            }

            Complex<N> omega = Omega(z);

            Complex<N> cos = Complex<N>.Cos(omega), sin = Complex<N>.Sin(omega);

            Complex<N> y = Complex<N>.Sqrt(2 / (MultiPrecision<N>.PI * z)) * (cos * c_even - sin * c_odd);

            return y;
        }

        public Complex<N> BesselY(Complex<N> z) {
            if (z.R.Sign == Sign.Minus) {
                return (cospi_nu, (z.I.Sign == Sign.Plus) ? -sinpi_nu : sinpi_nu) * BesselY(-z) +
                       (0, (z.I.Sign == Sign.Plus) ? 2 * cospi_nu : -2 * cospi_nu) * BesselJ(-z);
            }

            (Complex<N> c_even, Complex<N> c_odd) = BesselJYCoef(z);
            if (Complex<N>.IsNaN(c_even) || Complex<N>.IsNaN(c_odd)) {
                return Complex<N>.NaN;
            }

            Complex<N> omega = Omega(z);

            Complex<N> cos = Complex<N>.Cos(omega), sin = Complex<N>.Sin(omega);

            Complex<N> y = Complex<N>.Sqrt(2 / (MultiPrecision<N>.PI * z)) * (sin * c_even + cos * c_odd);

            return y;
        }

        public Complex<N> BesselI(Complex<N> z) {
            if (z.R.Sign == Sign.Minus) {
                return (cospi_nu, (z.I.Sign == Sign.Plus) ? sinpi_nu : -sinpi_nu) * BesselI(-z);
            }

            Complex<N> ci = BesselICoef(z), ck = BesselKCoef(z);
            if (Complex<N>.IsNaN(ci) || Complex<N>.IsNaN(ck)) {
                return Complex<N>.NaN;
            }

            Complex<N> y = Complex<N>.Sqrt(1 / (2 * MultiPrecision<N>.PI * z)) * (
                Complex<N>.Exp(z) * ci -
                (sinpi_nu, (z.I.Sign == Sign.Plus) ? -cospi_nu : cospi_nu) * Complex<N>.Exp(-z) * ck
            );

            return y;
        }

        public Complex<N> BesselK(Complex<N> z) {
            if (z.R.Sign == Sign.Minus) {
                return (cospi_nu, (z.I.Sign == Sign.Plus) ? -sinpi_nu : sinpi_nu) * BesselK(-z) -
                    (0, (z.I.Sign == Sign.Plus) ? MultiPrecision<N>.PI : -MultiPrecision<N>.PI) * BesselI(-z);
            }

            Complex<N> c = BesselKCoef(z);
            if (Complex<N>.IsNaN(c)) {
                return Complex<N>.NaN;
            }

            Complex<N> y = Complex<N>.Sqrt(MultiPrecision<N>.PI / (2 * z)) * Complex<N>.Exp(-z) * c;

            return y;
        }
    }
}
