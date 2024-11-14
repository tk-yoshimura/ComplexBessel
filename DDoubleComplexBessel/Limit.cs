using DoubleDouble;
using DoubleDoubleComplex;
using System.Collections.Concurrent;
using System.Diagnostics;

namespace DDoubleComplexBessel {
    public static class Limit {
        static readonly ConcurrentDictionary<ddouble, HankelExpansion> table = [];

        public static Complex BesselJ(ddouble nu, Complex z) {
            Debug.Assert(ddouble.IsPositive(z.R));
            Debug.Assert(ddouble.IsPositive(z.I));

            if (!table.TryGetValue(nu, out HankelExpansion hankel)) {
                hankel = new HankelExpansion(nu);
                table[nu] = hankel;
            }

            (Complex c_even, Complex c_odd) = hankel.BesselJYCoef(z);

            Complex omega = hankel.Omega(z);

            Complex cos = Complex.Cos(omega), sin = Complex.Sin(omega);

            Complex y = Complex.Sqrt(2d / (ddouble.Pi * z)) * (cos * c_even - sin * c_odd);

            return y;
        }

        public static Complex BesselY(ddouble nu, Complex z) {
            Debug.Assert(ddouble.IsPositive(z.R));
            Debug.Assert(ddouble.IsPositive(z.I));

            if (!table.TryGetValue(nu, out HankelExpansion hankel)) {
                hankel = new HankelExpansion(nu);
                table[nu] = hankel;
            }

            (Complex c_even, Complex c_odd) = hankel.BesselJYCoef(z);

            Complex omega = hankel.Omega(z);

            Complex cos = Complex.Cos(omega), sin = Complex.Sin(omega);

            Complex y = Complex.Sqrt(2d / (ddouble.Pi * z)) * (sin * c_even + cos * c_odd);

            return y;
        }

        public static Complex BesselI(ddouble nu, Complex z) {
            Debug.Assert(ddouble.IsPositive(z.R));
            Debug.Assert(ddouble.IsPositive(z.I));

            if (!table.TryGetValue(nu, out HankelExpansion hankel)) {
                hankel = new HankelExpansion(nu);
                table[nu] = hankel;
            }

            Complex ci = hankel.BesselICoef(z), ck = hankel.BesselKCoef(z);

            Complex y = Complex.Sqrt(1d / (2d * ddouble.Pi * z)) * (
                Complex.Exp(z) * ci -
                (SinCosPiCache.SinPi(nu), -SinCosPiCache.CosPi(nu)) * Complex.Exp(-z) * ck
            );

            return y;
        }

        public static Complex BesselK(ddouble nu, Complex z) {
            Debug.Assert(ddouble.IsPositive(nu));
            Debug.Assert(ddouble.IsPositive(z.R));
            Debug.Assert(ddouble.IsPositive(z.I));

            if (!table.TryGetValue(nu, out HankelExpansion hankel)) {
                hankel = new HankelExpansion(nu);
                table[nu] = hankel;
            }

            Complex c = hankel.BesselKCoef(z);

            Complex y = Complex.Sqrt(ddouble.Pi / (2d * z)) * Complex.Exp(-z) * c;

            return y;
        }

        private class HankelExpansion {
            public ddouble Nu { get; }

            private readonly List<ddouble> a_coef;

            public HankelExpansion(ddouble nu) {
                Nu = nu;
                a_coef = [1];
            }

            private ddouble ACoef(int n) {
                if (n < a_coef.Count) {
                    return a_coef[n];
                }

                lock (a_coef) {
                    for (int k = a_coef.Count; k <= n; k++) {
                        ddouble a = a_coef.Last() * (4d * Nu * Nu - checked((2 * k - 1) * (2 * k - 1))) / (k * 8);
                        a_coef.Add(a);
                    }

                    return a_coef[n];
                }
            }

            public Complex Omega(Complex z) {
                Complex omega = z - ddouble.Ldexp(2 * Nu + 1, -2) * ddouble.Pi;

                return omega;
            }

            public (Complex c_even, Complex c_odd) BesselJYCoef(Complex z, int max_term = 256) {
                Complex v = 1 / (z * z), w = -v;

                Complex c_even = ACoef(0), c_odd = ACoef(1);

                for (int k = 1; k <= max_term; k++) {
                    Complex dc_even = w * ACoef(2 * k);
                    Complex dc_odd = w * ACoef(2 * k + 1);

                    c_even += dc_even;
                    c_odd += dc_odd;

                    if (((long)Complex.ILogB(c_even) - Complex.ILogB(dc_even) >= 106L || Complex.IsZero(dc_even)) &&
                        ((long)Complex.ILogB(c_odd) - Complex.ILogB(dc_odd) >= 106L || Complex.IsZero(dc_odd))) {

                        IterationLogger.Log("BesselJY Hankel", k);

                        return (c_even, c_odd / z);
                    }

                    w *= -v;
                }

                return Complex.NaN;
            }

            public Complex BesselICoef(Complex z, int max_term = 256) {
                Complex v = 1d / z, w = -v;

                Complex c = ACoef(0);

                for (int k = 1; k <= max_term; k++) {
                    Complex dc = w * ACoef(k);

                    c += dc;

                    if ((long)Complex.ILogB(c) - Complex.ILogB(dc) >= 106L || Complex.IsZero(dc)) {
                        IterationLogger.Log("BesselI Hankel", k);

                        return c;
                    }

                    w *= -v;
                }

                return Complex.NaN;
            }

            public Complex BesselKCoef(Complex z, int max_term = 256) {
                Complex v = 1d / z, w = v;

                Complex c = ACoef(0);

                for (int k = 1; k <= max_term; k++) {
                    Complex dc = w * ACoef(k);

                    c += dc;

                    if ((long)Complex.ILogB(c) - Complex.ILogB(dc) >= 106L || Complex.IsZero(dc)) {
                        IterationLogger.Log("BesselK Hankel", k);

                        return c;
                    }

                    w *= v;
                }

                return Complex.NaN;
            }
        }
    }
}
