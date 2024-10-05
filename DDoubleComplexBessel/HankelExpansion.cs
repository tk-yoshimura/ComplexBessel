using DoubleDouble;
using DoubleDoubleComplex;
using System.Diagnostics;

namespace DDoubleComplexBessel {
    public class HankelExpansion {
        public ddouble Nu { get; }

        private readonly List<ddouble> a_coef;

        private readonly ddouble cospi_nu, sinpi_nu;

        public HankelExpansion(ddouble nu) {
            Nu = nu;
            a_coef = [1];
            cospi_nu = SinCosPICache.CosPI(Nu);
            sinpi_nu = SinCosPICache.SinPI(Nu);
        }

        public ddouble ACoef(int n) {
            for (int k = a_coef.Count; k <= n; k++) {
                ddouble a = a_coef.Last() * (4 * Nu * Nu - checked((2 * k - 1) * (2 * k - 1))) / (k * 8);
                a_coef.Add(a);
            }

            return a_coef[n];
        }

        private Complex Omega(Complex z) {
            Complex omega = z - ddouble.Ldexp(2 * Nu + 1, -2) * ddouble.PI;

            return omega;
        }

        private (Complex c_even, Complex c_odd) BesselJYCoef(Complex z, int max_term = 256) {
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

            return (ddouble.NaN, ddouble.NaN);
        }

        private Complex BesselICoef(Complex z, int max_term = 256) {
            Complex v = 1 / z, w = -v;

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

            return ddouble.NaN;
        }

        private Complex BesselKCoef(Complex z, int max_term = 256) {
            Complex v = 1 / z, w = v;

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

            return ddouble.NaN;
        }

        public Complex BesselJ(Complex z) {
            Debug.Assert(ddouble.IsPositive(z.R));
            Debug.Assert(ddouble.IsPositive(z.I));

            (Complex c_even, Complex c_odd) = BesselJYCoef(z);
            if (Complex.IsNaN(c_even) || Complex.IsNaN(c_odd)) {
                return Complex.NaN;
            }

            Complex omega = Omega(z);

            Complex cos = Complex.Cos(omega), sin = Complex.Sin(omega);

            Complex y = Complex.Sqrt(2 / (ddouble.PI * z)) * (cos * c_even - sin * c_odd);

            return y;
        }

        public Complex BesselY(Complex z) {
            Debug.Assert(ddouble.IsPositive(z.R));
            Debug.Assert(ddouble.IsPositive(z.I));

            (Complex c_even, Complex c_odd) = BesselJYCoef(z);
            if (Complex.IsNaN(c_even) || Complex.IsNaN(c_odd)) {
                return Complex.NaN;
            }

            Complex omega = Omega(z);

            Complex cos = Complex.Cos(omega), sin = Complex.Sin(omega);

            Complex y = Complex.Sqrt(2 / (ddouble.PI * z)) * (sin * c_even + cos * c_odd);

            return y;
        }

        public Complex BesselI(Complex z) {
            Debug.Assert(ddouble.IsPositive(z.R));
            Debug.Assert(ddouble.IsPositive(z.I));

            Complex ci = BesselICoef(z), ck = BesselKCoef(z);
            if (Complex.IsNaN(ci) || Complex.IsNaN(ck)) {
                return Complex.NaN;
            }

            Complex y = Complex.Sqrt(1 / (2 * ddouble.PI * z)) * (
                Complex.Exp(z) * ci -
                (sinpi_nu, -cospi_nu) * Complex.Exp(-z) * ck
            );

            return y;
        }

        public Complex BesselK(Complex z) {
            Debug.Assert(ddouble.IsPositive(Nu));
            Debug.Assert(ddouble.IsPositive(z.R));
            Debug.Assert(ddouble.IsPositive(z.I));

            Complex c = BesselKCoef(z);
            if (Complex.IsNaN(c)) {
                return Complex.NaN;
            }

            Complex y = Complex.Sqrt(ddouble.PI / (2 * z)) * Complex.Exp(-z) * c;

            return y;
        }
    }
}
