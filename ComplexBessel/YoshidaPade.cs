using MultiPrecision;
using MultiPrecisionComplex;
using System.Collections.ObjectModel;

namespace ComplexBessel {

    public static class YoshidaPade<N> where N : struct, IConstant {
        const int m = 56;

        private static readonly ReadOnlyCollection<ReadOnlyCollection<MultiPrecision<N>>> ess_coef_table;
        private static readonly Dictionary<MultiPrecision<N>, ReadOnlyCollection<(MultiPrecision<N> c, MultiPrecision<N> s)>> cds_coef_table = [];

        static YoshidaPade() {
            MultiPrecision<N>[][] dss_mp4 = YoshidaCoef<N>.Table(m);

            MultiPrecision<N> n = dss_mp4[0][0];

            ReadOnlyCollection<MultiPrecision<N>>[] dss = new ReadOnlyCollection<MultiPrecision<N>>[dss_mp4.Length];

            for (int i = 0; i < dss_mp4.Length; i++) {
                dss[i] = new(dss_mp4[i].Select(v => v / n).ToArray());
            }

            ess_coef_table = new(dss);
        }

        public static Complex<N> BesselK(MultiPrecision<N> nu, Complex<N> z) {
            if (nu < 2d) {
                if (!cds_coef_table.TryGetValue(nu, out ReadOnlyCollection<(MultiPrecision<N> c, MultiPrecision<N> s)> cds_table)) {
                    cds_table = Table(nu);
                    cds_coef_table.Add(nu, cds_table);
                }

                ReadOnlyCollection<(MultiPrecision<N>, MultiPrecision<N>)> cds = cds_table;

                Complex<N> y = Value(z, cds);

                return y;
            }
            else {
                int n = (int)MultiPrecision<N>.Floor(nu);
                MultiPrecision<N> alpha = nu - n;

                Complex<N> y0 = BesselK(alpha, z);
                Complex<N> y1 = BesselK(alpha + 1d, z);

                Complex<N> v = 1d / z;

                for (int k = 1; k < n; k++) {
                    (y1, y0) = (Complex<N>.Ldexp(k + alpha, 1) * v * y1 + y0, y1);
                }

                return y1;
            }
        }

        private static Complex<N> Value(Complex<N> z, ReadOnlyCollection<(MultiPrecision<N> c, MultiPrecision<N> d)> cds) {
            Complex<N> t = 1d / z;
            (Complex<N> sc, Complex<N> sd) = cds[0];

            for (int i = 1; i < cds.Count; i++) {
                (MultiPrecision<N> c, MultiPrecision<N> d) = cds[i];

                sc = sc * t + c;
                sd = sd * t + d;
            }

            Complex<N> y = Complex<N>.Sqrt(Complex<N>.Ldexp(t * MultiPrecision<N>.PI, -1)) * sc / sd;

            y *= Complex<N>.Exp(-z);

            return y;
        }

        private static ReadOnlyCollection<(MultiPrecision<N> c, MultiPrecision<N> d)> Table(MultiPrecision<N> nu) {
            int m = ess_coef_table.Count - 1;

            MultiPrecision<N> squa_nu = nu * nu;
            List<(MultiPrecision<N> c, MultiPrecision<N> d)> cds = [];
            MultiPrecision<N>[] us = new MultiPrecision<N>[m + 1], vs = new MultiPrecision<N>[m];

            MultiPrecision<N> u = 1d;
            for (int i = 0; i <= m; i++) {
                us[i] = u;
                u *= squa_nu;
            }
            for (int i = 0; i < m; i++) {
                MultiPrecision<N> r = m - i + 0.5d;
                vs[i] = r * r - squa_nu;
            }

            for (int i = 0; i <= m; i++) {
                ReadOnlyCollection<MultiPrecision<N>> es = ess_coef_table[i];
                MultiPrecision<N> d = es[i], c = 0d;

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
