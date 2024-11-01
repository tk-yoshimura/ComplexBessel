using MultiPrecision;
using MultiPrecisionComplex;
using System.Diagnostics;

namespace ComplexBessel {
    public static class Limit<N> where N : struct, IConstant {
        static readonly Dictionary<MultiPrecision<N>, HankelExpansion> table = [];

        public static Complex<N> BesselJ(MultiPrecision<N> nu, Complex<N> z) {
            Debug.Assert(MultiPrecision<N>.IsPositive(z.R));
            Debug.Assert(MultiPrecision<N>.IsPositive(z.I));

            if (!table.TryGetValue(nu, out HankelExpansion hankel)) {
                hankel = new HankelExpansion(nu);
                table.Add(nu, hankel);
            }

            (Complex<N> c_even, Complex<N> c_odd) = hankel.BesselJYCoef(z);

            Complex<N> omega = hankel.Omega(z);

            Complex<N> cos = Complex<N>.Cos(omega), sin = Complex<N>.Sin(omega);

            Complex<N> y = Complex<N>.Sqrt(2d / (MultiPrecision<N>.Pi * z)) * (cos * c_even - sin * c_odd);

            return y;
        }

        public static Complex<N> BesselY(MultiPrecision<N> nu, Complex<N> z) {
            Debug.Assert(MultiPrecision<N>.IsPositive(z.R));
            Debug.Assert(MultiPrecision<N>.IsPositive(z.I));

            if (!table.TryGetValue(nu, out HankelExpansion hankel)) {
                hankel = new HankelExpansion(nu);
                table.Add(nu, hankel);
            }


            (Complex<N> c_even, Complex<N> c_odd) = hankel.BesselJYCoef(z);

            Complex<N> omega = hankel.Omega(z);

            Complex<N> cos = Complex<N>.Cos(omega), sin = Complex<N>.Sin(omega);

            Complex<N> y = Complex<N>.Sqrt(2d / (MultiPrecision<N>.Pi * z)) * (sin * c_even + cos * c_odd);

            return y;
        }

        public static Complex<N> BesselI(MultiPrecision<N> nu, Complex<N> z) {
            Debug.Assert(MultiPrecision<N>.IsPositive(z.R));
            Debug.Assert(MultiPrecision<N>.IsPositive(z.I));

            if (!table.TryGetValue(nu, out HankelExpansion hankel)) {
                hankel = new HankelExpansion(nu);
                table.Add(nu, hankel);
            }

            Complex<N> ci = hankel.BesselICoef(z), ck = hankel.BesselKCoef(z);

            Complex<N> y = Complex<N>.Sqrt(1d / (2d * MultiPrecision<N>.Pi * z)) * (
                Complex<N>.Exp(z) * ci -
                (SinCosPiCache<N>.SinPi(nu), -SinCosPiCache<N>.CosPi(nu)) * Complex<N>.Exp(-z) * ck
            );

            return y;
        }

        public static Complex<N> BesselK(MultiPrecision<N> nu, Complex<N> z) {
            Debug.Assert(MultiPrecision<N>.IsPositive(nu));
            Debug.Assert(MultiPrecision<N>.IsPositive(z.R));
            Debug.Assert(MultiPrecision<N>.IsPositive(z.I));

            if (!table.TryGetValue(nu, out HankelExpansion hankel)) {
                hankel = new HankelExpansion(nu);
                table.Add(nu, hankel);
            }

            Complex<N> c = hankel.BesselKCoef(z);

            Complex<N> y = Complex<N>.Sqrt(MultiPrecision<N>.Pi / (2d * z)) * Complex<N>.Exp(-z) * c;

            return y;
        }

        public class HankelExpansion {
            public MultiPrecision<N> Nu { get; }

            private readonly List<MultiPrecision<N>> a_coef;

            public HankelExpansion(MultiPrecision<N> nu) {
                Nu = nu;
                a_coef = [1];
            }

            public MultiPrecision<N> ACoef(int n) {
                for (int k = a_coef.Count; k <= n; k++) {
                    MultiPrecision<N> a = a_coef.Last() * (4d * Nu * Nu - checked((2 * k - 1) * (2 * k - 1))) / (k * 8);
                    a_coef.Add(a);
                }

                return a_coef[n];
            }

            public Complex<N> Omega(Complex<N> z) {
                Complex<N> omega = z - MultiPrecision<N>.Ldexp(2 * Nu + 1, -2) * MultiPrecision<N>.Pi;

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

                    if ((c_even.Exponent - dc_even.Exponent >= MultiPrecision<N>.Bits || Complex<N>.IsZero(dc_even)) &&
                        (c_odd.Exponent - dc_odd.Exponent >= MultiPrecision<N>.Bits || Complex<N>.IsZero(dc_odd))) {

                        return (c_even, c_odd / z);
                    }

                    w *= -v;
                }

                return Complex<N>.NaN;
            }

            public Complex<N> BesselICoef(Complex<N> z, int max_term = 256) {
                Complex<N> v = 1d / z, w = -v;

                Complex<N> c = ACoef(0);

                for (int k = 1; k <= max_term; k++) {
                    Complex<N> dc = w * ACoef(k);

                    c += dc;

                    if (c.Exponent - dc.Exponent >= MultiPrecision<N>.Bits || Complex<N>.IsZero(dc)) {
                        return c;
                    }

                    w *= -v;
                }

                return Complex<N>.NaN;
            }

            public Complex<N> BesselKCoef(Complex<N> z, int max_term = 256) {
                Complex<N> v = 1d / z, w = v;

                Complex<N> c = ACoef(0);

                for (int k = 1; k <= max_term; k++) {
                    Complex<N> dc = w * ACoef(k);

                    c += dc;

                    if (c.Exponent - dc.Exponent >= MultiPrecision<N>.Bits || Complex<N>.IsZero(dc)) {
                        return c;
                    }

                    w *= v;
                }

                return Complex<N>.NaN;
            }
        }
    }
}
