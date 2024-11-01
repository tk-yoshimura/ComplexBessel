using ComplexBessel;
using DoubleDouble;
using DoubleDoubleComplex;
using MultiPrecision;
using System.Collections.ObjectModel;

namespace DDoubleComplexBessel {

    public static class YoshidaPade {
        const int m = 38;

        private static readonly ReadOnlyCollection<ReadOnlyCollection<ddouble>> ess_coef_table;
        private static readonly Dictionary<ddouble, ReadOnlyCollection<(ddouble c, ddouble s)>> cds_coef_table = [];

        static YoshidaPade() {
            MultiPrecision<Pow2.N4>[][] dss_mp4 = YoshidaCoef<Pow2.N4>.Table(m);

            MultiPrecision<Pow2.N4> n = dss_mp4[0][0];

            ReadOnlyCollection<ddouble>[] dss = new ReadOnlyCollection<ddouble>[dss_mp4.Length];

            for (int i = 0; i < dss_mp4.Length; i++) {
                dss[i] = new(dss_mp4[i].Select(v => (ddouble)(v / n).ToString()).ToArray());
            }

            ess_coef_table = new(dss);
        }

        public static Complex BesselK(ddouble nu, Complex z) {
            if (nu < 2d) {
                if (!cds_coef_table.TryGetValue(nu, out ReadOnlyCollection<(ddouble c, ddouble s)> cds_table)) {
                    cds_table = Table(nu);
                    cds_coef_table.Add(nu, cds_table);
                }

                ReadOnlyCollection<(ddouble, ddouble)> cds = cds_table;

                Complex y = Value(z, cds);

                return y;
            }
            else {
                int n = (int)ddouble.Floor(nu);
                ddouble alpha = nu - n;

                Complex y0 = BesselK(alpha, z);
                Complex y1 = BesselK(alpha + 1d, z);

                Complex v = 1d / z;

                for (int k = 1; k < n; k++) {
                    (y1, y0) = (Complex.Ldexp(k + alpha, 1) * v * y1 + y0, y1);
                }

                return y1;
            }
        }

        private static Complex Value(Complex z, ReadOnlyCollection<(ddouble c, ddouble d)> cds) {
            Complex t = 1d / z;
            (Complex sc, Complex sd) = cds[0];

            for (int i = 1; i < cds.Count; i++) {
                (ddouble c, ddouble d) = cds[i];

                sc = sc * t + c;
                sd = sd * t + d;
            }

            Complex y = Complex.Sqrt(Complex.Ldexp(t * ddouble.Pi, -1)) * sc / sd;

            y *= Complex.Exp(-z);

            return y;
        }

        private static ReadOnlyCollection<(ddouble c, ddouble d)> Table(ddouble nu) {
            int m = ess_coef_table.Count - 1;

            ddouble squa_nu = nu * nu;
            List<(ddouble c, ddouble d)> cds = [];
            ddouble[] us = new ddouble[m + 1], vs = new ddouble[m];

            ddouble u = 1d;
            for (int i = 0; i <= m; i++) {
                us[i] = u;
                u *= squa_nu;
            }
            for (int i = 0; i < m; i++) {
                ddouble r = m - i + 0.5d;
                vs[i] = r * r - squa_nu;
            }

            for (int i = 0; i <= m; i++) {
                ReadOnlyCollection<ddouble> es = ess_coef_table[i];
                ddouble d = es[i], c = 0d;

                for (int l = 0; l < i; l++) {
                    d *= vs[l];
                }
                for (int j = 0; j <= i; j++) {
                    c += es[j] * us[j];
                }

                cds.Add((c, d));
            }

            cds.Reverse();

            return Array.AsReadOnly(cds.ToArray());
        }
    }
}
