using ComplexBessel;
using DoubleDouble;
using DoubleDoubleComplex;
using MultiPrecision;
using System.Diagnostics;

namespace DDoubleComplexBessel {
    public class YoshidaPade {
        const int m = 36;
        private static readonly ddouble[][] dss;

        public ddouble Nu { get; }
        private readonly ddouble[] cs0, ds0, cs1, ds1;

        static YoshidaPade() {
            MultiPrecision<Pow2.N4>[][] dss_mp4 = YoshidaCoef<Pow2.N4>.Table(m);

            MultiPrecision<Pow2.N4> n = dss_mp4[0][0];

            dss = new ddouble[dss_mp4.Length][];

            for (int i = 0; i < dss.Length; i++) {
                dss[i] = dss_mp4[i].Select(v => (ddouble)(v / n).ToString()).ToArray();
            }
        }

        public YoshidaPade(ddouble nu) {
            Nu = ddouble.Abs(nu);

            (cs0, ds0) = Table(Nu - ddouble.Floor(Nu), dss);
            (cs1, ds1) = Table(Nu - ddouble.Floor(Nu) + 1, dss);
        }

        public static (ddouble[] cs, ddouble[] ds) Table(ddouble nu, ddouble[][] dss) {
            int m = dss.Length - 1;

            ddouble squa_nu = nu * nu;
            ddouble[] cs = new ddouble[m + 1], ds = new ddouble[m + 1];

            for (int i = 0; i <= m; i++) {
                ddouble d = dss[i][i], c = 0;

                for (int l = 0; l < i; l++) {
                    d *= ddouble.Square(m - l + ddouble.Point5) - squa_nu;
                }

                ddouble u = 1;
                for (int j = 0; j <= i; j++) {
                    c += dss[i][j] * u;
                    u *= squa_nu;
                }

                cs[i] = c;
                ds[i] = d;
            }

            return (cs, ds);
        }

        public Complex BesselK(Complex z) {
            Debug.Assert(ddouble.IsPositive(z.R));
            Debug.Assert(ddouble.IsPositive(z.I));

            if (Nu < 1d) {
                Complex y = Value(z, cs0, ds0);

                return y;
            }
            else if (Nu < 2d) {
                Complex y = Value(z, cs1, ds1);

                return y;
            }
            else {
                int n = (int)ddouble.Floor(Nu);
                ddouble alpha = Nu - n;

                Complex y0 = Value(z, cs0, ds0);
                Complex y1 = Value(z, cs1, ds1);

                Complex v = 1d / z;

                for (int k = 1; k < n; k++) {
                    (y1, y0) = (Complex.Ldexp(k + alpha, 1) * v * y1 + y0, y1);
                }

                return y1;
            }
        }

        static Complex Value(Complex z, ddouble[] cs, ddouble[] ds) {
            Complex t = 1 / z, tn = 1;
            Complex c = 0, d = 0;

            for (int j = 0; j < cs.Length; j++) {
                c += cs[j] * tn;
                d += ds[j] * tn;
                tn *= t;
            }

            Complex y = Complex.Sqrt(t * (ddouble.PI / 2)) * c / d;
            y *= Complex.Exp(-z);

            return y;
        }
    }
}
