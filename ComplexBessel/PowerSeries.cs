using MultiPrecision;
using MultiPrecisionComplex;
using System.Collections.ObjectModel;
using System.Diagnostics;

namespace ComplexBessel {
    public class PowerSeries<N> where N : struct, IConstant {
        private static readonly Dictionary<MultiPrecision<N>, DoubleFactDenomTable> dfactdenom_coef_table = [];
        private static readonly Dictionary<MultiPrecision<N>, X2DenomTable> x2denom_coef_table = [];
        private static readonly Dictionary<MultiPrecision<N>, GammaDenomTable> gammadenom_coef_table = [];
        private static readonly Dictionary<MultiPrecision<N>, GammaTable> gamma_coef_table = [];
        private static readonly Dictionary<MultiPrecision<N>, GammaPNTable> gammapn_coef_table = [];
        private static readonly YCoefTable y_coef_table = new();
        private static readonly Y0CoefTable y0_coef_table = new();
        private static readonly Dictionary<int, YNCoefTable> yn_coef_table = [];
        private static readonly Dictionary<int, ReadOnlyCollection<MultiPrecision<N>>> yn_finitecoef_table = [];
        private static readonly KCoefTable k_coef_table = new();
        private static readonly K0CoefTable k0_coef_table = new();
        private static readonly K1CoefTable k1_coef_table = new();

        public static Complex<N> BesselJ(MultiPrecision<N> nu, Complex<N> z) {
            Debug.Assert(MultiPrecision<N>.IsPositive(z.R));
            Debug.Assert(MultiPrecision<N>.IsPositive(z.I));

            if (MultiPrecision<N>.IsNegative(nu) && BesselUtil<N>.NearlyInteger(nu, out int n)) {
                Complex<N> y = BesselJ(-nu, z);

                return (n & 1) == 0 ? y : -y;
            }
            else {
                Complex<N> y = BesselJKernel(nu, z, terms: 256);

                return y;
            }
        }

        public static Complex<N> BesselY(MultiPrecision<N> nu, Complex<N> z) {
            Debug.Assert(MultiPrecision<N>.IsPositive(z.R));
            Debug.Assert(MultiPrecision<N>.IsPositive(z.I));

            if (BesselUtil<N>.NearlyInteger(nu, out int n)) {
                Complex<N> y = BesselYKernel(n, z, terms: 256);

                return y;
            }
            else {
                Complex<N> y = BesselYKernel(nu, z, terms: 256);

                return y;
            }
        }

        public static Complex<N> BesselI(MultiPrecision<N> nu, Complex<N> z) {
            Debug.Assert(MultiPrecision<N>.IsPositive(z.R));
            Debug.Assert(MultiPrecision<N>.IsPositive(z.I));

            if (MultiPrecision<N>.IsNegative(nu) && BesselUtil<N>.NearlyInteger(nu, out _)) {
                Complex<N> y = BesselI(-nu, z);

                return y;
            }
            else {
                Complex<N> y = BesselIKernel(nu, z, terms: 256);

                return y;
            }
        }

        public static Complex<N> BesselK(MultiPrecision<N> nu, Complex<N> z) {
            Debug.Assert(MultiPrecision<N>.IsPositive(z.R));
            Debug.Assert(MultiPrecision<N>.IsPositive(z.I));

            if (BesselUtil<N>.NearlyInteger(nu, out int n)) {
                Complex<N> y = BesselKKernel(n, z, terms: 256);

                return y;
            }
            else {
                Complex<N> y = BesselKKernel(nu, z, terms: 256);

                return y;
            }
        }

        private static Complex<N> BesselJKernel(MultiPrecision<N> nu, Complex<N> z, int terms) {
            if (!dfactdenom_coef_table.TryGetValue(nu, out DoubleFactDenomTable r)) {
                r = new DoubleFactDenomTable(nu);
                dfactdenom_coef_table.Add(nu, r);
            }
            if (!x2denom_coef_table.TryGetValue(nu, out X2DenomTable d)) {
                d = new X2DenomTable(nu);
                x2denom_coef_table.Add(nu, d);
            }

            Complex<N> z2 = z * z, z4 = z2 * z2;

            Complex<N> c = 0d, u = Complex<N>.Pow(Complex<N>.Ldexp(z, -1), nu);

            for (int k = 0; k <= terms; k++) {
                c = SeriesUtil<N>.Add(c, u * r[k], 1d, -z2 * d[k], out bool convergence);

                if (convergence) {
                    break;
                }

                u *= z4;

                if (!Complex<N>.IsFinite(c)) {
                    break;
                }
            }

            return c;
        }

        private static Complex<N> BesselYKernel(MultiPrecision<N> nu, Complex<N> z, int terms) {
            if (!gamma_coef_table.TryGetValue(nu, out GammaTable g)) {
                g = new GammaTable(nu);
                gamma_coef_table.Add(nu, g);
            }
            if (!gammapn_coef_table.TryGetValue(nu, out GammaPNTable gpn)) {
                gpn = new GammaPNTable(nu);
                gammapn_coef_table.Add(nu, gpn);
            }

            YCoefTable r = y_coef_table;

            MultiPrecision<N> cos = SinCosPICache<N>.CosPI(nu), sin = SinCosPICache<N>.SinPI(nu);
            Complex<N> p = Complex<N>.IsZero(cos) ? 0d : Complex<N>.Pow(z, Complex<N>.Ldexp(nu, 1)) * cos, s = Complex<N>.Ldexp(Complex<N>.Pow(Complex<N>.Ldexp(z, 1), nu), 2);

            Complex<N> z2 = z * z, z4 = z2 * z2;

            Complex<N> c = 0d, u = 1d / sin;

            for (int k = 0, t = 1; k <= terms; k++, t += 2) {
                Complex<N> a = t * s * g[t], q = gpn[t];
                Complex<N> pa = p / a, qa = q / a;

                Complex<N> v = 4 * t * t - z2;
                c = SeriesUtil<N>.Add(c, u * r[k], 4 * t * nu * (pa + qa), v * (pa - qa), out bool convergence);

                if (convergence && v.Exponent >= -4) {
                    break;
                }

                u *= z4;

                if (!Complex<N>.IsFinite(c)) {
                    break;
                }
            }

            return c;
        }

        private static Complex<N> BesselYKernel(int n, Complex<N> z, int terms) {
            if (n < 0) {
                Complex<N> y = BesselYKernel(-n, z, terms);

                return (n & 1) == 0 ? y : -y;
            }
            else if (n == 0) {
                return BesselY0Kernel(z, terms);
            }
            else if (n == 1) {
                return BesselY1Kernel(z, terms);
            }
            else {
                return BesselYNKernel(n, z, terms);
            }
        }

        private static Complex<N> BesselY0Kernel(Complex<N> z, int terms) {
            if (!dfactdenom_coef_table.TryGetValue(0, out DoubleFactDenomTable r)) {
                r = new DoubleFactDenomTable(0);
                dfactdenom_coef_table.Add(0, r);
            }
            if (!x2denom_coef_table.TryGetValue(0, out X2DenomTable d)) {
                d = new X2DenomTable(0);
                x2denom_coef_table.Add(0, d);
            }

            Y0CoefTable q = y0_coef_table;

            Complex<N> h = Complex<N>.Log(Complex<N>.Ldexp(z, -1)) + MultiPrecision<N>.EulerGamma;

            Complex<N> z2 = z * z, z4 = z2 * z2;

            Complex<N> c = 0d, u = Complex<N>.Ldexp(MultiPrecision<N>.RcpPI, 1);

            for (int k = 0; k <= terms; k++) {
                Complex<N> s = u * r[k], t = h - MultiPrecision<N>.HarmonicNumber(2 * k);
                c = SeriesUtil<N>.Add(c, s * t, 1d, -z2 * d[k], out bool convergence1);
                c = SeriesUtil<N>.Add(c, s, z2 * q[k], out bool convergence2);

                if (convergence1 && convergence2 && t.Exponent >= -4) {
                    break;
                }

                u *= z4;
            }

            return c;
        }

        private static Complex<N> BesselY1Kernel(Complex<N> z, int terms) {
            if (!dfactdenom_coef_table.TryGetValue(1, out DoubleFactDenomTable r)) {
                r = new DoubleFactDenomTable(1);
                dfactdenom_coef_table.Add(1, r);
            }
            if (!x2denom_coef_table.TryGetValue(1, out X2DenomTable d)) {
                d = new X2DenomTable(1);
                x2denom_coef_table.Add(1, d);
            }
            if (!yn_coef_table.TryGetValue(1, out YNCoefTable q)) {
                q = new YNCoefTable(1);
                yn_coef_table.Add(1, q);
            }

            Complex<N> h = Complex<N>.Ldexp(Complex<N>.Log(Complex<N>.Ldexp(z, -1)) + MultiPrecision<N>.EulerGamma, 1);

            Complex<N> z2 = z * z, z4 = z2 * z2;

            Complex<N> c = -2d / (z * MultiPrecision<N>.PI), u = z / Complex<N>.Ldexp(MultiPrecision<N>.PI, 1);

            for (int k = 0; k <= terms; k++) {
                Complex<N> s = u * r[k], t = h - MultiPrecision<N>.HarmonicNumber(2 * k) - MultiPrecision<N>.HarmonicNumber(2 * k + 1);
                c = SeriesUtil<N>.Add(c, s * t, 1d, -z2 * d[k], out bool convergence1);
                c = SeriesUtil<N>.Add(c, s, z2 * q[k], out bool convergence2);

                if (convergence1 && convergence2 && t.Exponent >= -4) {
                    break;
                }

                u *= z4;
            }

            return c;
        }

        private static Complex<N> BesselYNKernel(int n, Complex<N> z, int terms) {
            if (!dfactdenom_coef_table.TryGetValue(n, out DoubleFactDenomTable r)) {
                r = new DoubleFactDenomTable(n);
                dfactdenom_coef_table.Add(n, r);
            }
            if (!x2denom_coef_table.TryGetValue(n, out X2DenomTable d)) {
                d = new X2DenomTable(n);
                x2denom_coef_table.Add(n, d);
            }
            if (!yn_coef_table.TryGetValue(n, out YNCoefTable q)) {
                q = new YNCoefTable(n);
                yn_coef_table.Add(n, q);
            }
            if (!yn_finitecoef_table.TryGetValue(n, out ReadOnlyCollection<MultiPrecision<N>> f)) {
                f = YNFiniteCoefTable.Value(n);
                yn_finitecoef_table.Add(n, f);
            }

            Complex<N> c = 0d;
            Complex<N> z2 = z * z, z4 = z2 * z2;
            Complex<N> u = 1d, v = 1d, w = Complex<N>.Ldexp(z2, -2);

            for (int k = 0; k < n; k++) {
                c += v * f[k];
                v *= w;
            }
            c /= -v;

            Complex<N> h = Complex<N>.Ldexp(Complex<N>.Log(Complex<N>.Ldexp(z, -1)) + MultiPrecision<N>.EulerGamma, 1);

            for (int k = 0; k <= terms; k++) {
                Complex<N> s = u * r[k], t = (h - MultiPrecision<N>.HarmonicNumber(2 * k) - MultiPrecision<N>.HarmonicNumber(2 * k + n));
                c = SeriesUtil<N>.Add(c, s * t, 1d, -z2 * d[k], out bool convergence1);
                c = SeriesUtil<N>.Add(c, s, z2 * q[k], out bool convergence2);

                if (convergence1 && convergence2 && t.Exponent >= -4) {
                    break;
                }

                u *= z4;
            }

            Complex<N> y = c * MultiPrecision<N>.RcpPI * Complex<N>.Pow(Complex<N>.Ldexp(z, -1), n);

            return y;
        }

        private static Complex<N> BesselIKernel(MultiPrecision<N> nu, Complex<N> z, int terms) {
            if (!dfactdenom_coef_table.TryGetValue(nu, out DoubleFactDenomTable r)) {
                r = new DoubleFactDenomTable(nu);
                dfactdenom_coef_table.Add(nu, r);
            }
            if (!x2denom_coef_table.TryGetValue(nu, out X2DenomTable d)) {
                d = new X2DenomTable(nu);
                x2denom_coef_table.Add(nu, d);
            }

            Complex<N> z2 = z * z, z4 = z2 * z2;

            Complex<N> c = 0d, u = Complex<N>.Pow(Complex<N>.Ldexp(z, -1), nu);

            for (int k = 0; k <= terms; k++) {
                c = SeriesUtil<N>.Add(c, u * r[k], 1d, z2 * d[k], out bool convergence);

                if (convergence) {
                    break;
                }

                u *= z4;

                if (!Complex<N>.IsFinite(c)) {
                    break;
                }
            }

            return c;
        }

        private static Complex<N> BesselKKernel(MultiPrecision<N> nu, Complex<N> z, int terms) {
            if (!gammadenom_coef_table.TryGetValue(nu, out GammaDenomTable gp)) {
                gp = new GammaDenomTable(nu);
                gammadenom_coef_table.Add(nu, gp);
            }
            if (!gammadenom_coef_table.TryGetValue(-nu, out GammaDenomTable gn)) {
                gn = new GammaDenomTable(-nu);
                gammadenom_coef_table.Add(-nu, gn);
            }

            KCoefTable r = k_coef_table;

            Complex<N> tp = Complex<N>.Pow(Complex<N>.Ldexp(z, -1), nu), tn = 1d / tp;

            Complex<N> z2 = z * z;

            Complex<N> c = 0d, u = MultiPrecision<N>.PI / Complex<N>.Ldexp(MultiPrecision<N>.SinPI(nu), 1);

            for (int k = 0; k <= terms; k++) {
                c = SeriesUtil<N>.Add(c, u * r[k], tn * gn[k], -tp * gp[k], out bool convergence);

                if (convergence) {
                    break;
                }

                u *= z2;

                if (!Complex<N>.IsFinite(c)) {
                    break;
                }
            }

            return c;
        }

        private static Complex<N> BesselKKernel(int n, Complex<N> z, int terms) {
            if (n == 0) {
                return BesselK0Kernel(z, terms);
            }
            else if (n == 1) {
                return BesselK1Kernel(z, terms);
            }
            else {
                return BesselKNKernel(n, z, terms);
            }
        }

        private static Complex<N> BesselK0Kernel(Complex<N> z, int terms) {
            K0CoefTable r = k0_coef_table;
            Complex<N> h = -Complex<N>.Log(Complex<N>.Ldexp(z, -1)) - MultiPrecision<N>.EulerGamma;

            Complex<N> z2 = z * z;

            Complex<N> c = 0d, u = 1d;

            for (int k = 0; k <= terms; k++) {
                c = SeriesUtil<N>.Add(c, u * r[k], h, MultiPrecision<N>.HarmonicNumber(k), out bool convergence);

                if (convergence) {
                    break;
                }

                u *= z2;
            }

            return c;
        }

        private static Complex<N> BesselK1Kernel(Complex<N> z, int terms) {
            K1CoefTable r = k1_coef_table;
            Complex<N> h = Complex<N>.Log(Complex<N>.Ldexp(z, -1)) + MultiPrecision<N>.EulerGamma;

            Complex<N> z2 = z * z;

            Complex<N> c = 1d / z, u = Complex<N>.Ldexp(z, -1);

            for (int k = 0; k <= terms; k++) {
                c = SeriesUtil<N>.Add(c, u * r[k], h, -MultiPrecision<N>.Ldexp(MultiPrecision<N>.HarmonicNumber(k) + MultiPrecision<N>.HarmonicNumber(k + 1), -1), out bool convergence);

                if (convergence) {
                    break;
                }

                u *= z2;
            }

            return c;
        }

        private static Complex<N> BesselKNKernel(int n, Complex<N> z, int terms) {
            Complex<N> v = 1d / z;
            Complex<N> y0 = BesselK0Kernel(z, terms);
            Complex<N> y1 = BesselK1Kernel(z, terms);

            for (int k = 1; k < n; k++) {
                (y1, y0) = (2 * k * v * y1 + y0, y1);
            }

            return y1;
        }

        private class DoubleFactDenomTable {
            private MultiPrecision<N> c;
            private readonly MultiPrecision<N> nu;
            private readonly List<MultiPrecision<N>> table = [];

            public DoubleFactDenomTable(MultiPrecision<N> nu) {
                this.c = MultiPrecision<N>.Gamma(nu + 1d);
                this.nu = nu;
                this.table.Add(MultiPrecision<N>.Rcp(c));
            }

            public MultiPrecision<N> this[int k] => Value(k);

            public MultiPrecision<N> Value(int k) {
                Debug.Assert(k >= 0);

                if (k < table.Count) {
                    return table[k];
                }

                for (long i = table.Count; i <= k; i++) {
                    c *= checked((nu + 2 * i) * (nu + (2 * i - 1)) * (32 * i * (2 * i - 1)));

                    table.Add(MultiPrecision<N>.Rcp(c));
                }

                return table[k];
            }
        }

        private class X2DenomTable {
            private readonly MultiPrecision<N> nu;
            private readonly List<MultiPrecision<N>> table = [];

            public X2DenomTable(MultiPrecision<N> nu) {
                MultiPrecision<N> a = MultiPrecision<N>.Rcp(4d * (nu + 1d));

                this.nu = nu;
                this.table.Add(a);
            }

            public MultiPrecision<N> this[int k] => Value(k);

            public MultiPrecision<N> Value(int k) {
                Debug.Assert(k >= 0);

                if (k < table.Count) {
                    return table[k];
                }

                for (long i = table.Count; i <= k; i++) {
                    MultiPrecision<N> a = MultiPrecision<N>.Rcp(checked(4d * (2 * i + 1) * (2 * i + 1 + nu)));

                    table.Add(a);
                }

                return table[k];
            }
        }

        private class GammaDenomTable {
            private MultiPrecision<N> c;
            private readonly MultiPrecision<N> nu;
            private readonly List<MultiPrecision<N>> table = [];

            public GammaDenomTable(MultiPrecision<N> nu) {
                this.c = MultiPrecision<N>.Gamma(nu + 1d);
                this.nu = nu;
                this.table.Add(MultiPrecision<N>.Rcp(c));
            }

            public MultiPrecision<N> this[int k] => Value(k);

            public MultiPrecision<N> Value(int k) {
                Debug.Assert(k >= 0);

                if (k < table.Count) {
                    return table[k];
                }

                for (int i = table.Count; i <= k; i++) {
                    c *= nu + i;

                    table.Add(MultiPrecision<N>.Rcp(c));
                }

                return table[k];
            }
        }

        private class GammaTable {
            private MultiPrecision<N> c;
            private readonly MultiPrecision<N> nu;
            private readonly List<MultiPrecision<N>> table = [];

            public GammaTable(MultiPrecision<N> nu) {
                this.c = MultiPrecision<N>.Gamma(nu + 1d);
                this.nu = nu;
                this.table.Add(c);
            }

            public MultiPrecision<N> this[int k] => Value(k);

            public MultiPrecision<N> Value(int k) {
                Debug.Assert(k >= 0);

                if (k < table.Count) {
                    return table[k];
                }

                for (int i = table.Count; i <= k; i++) {
                    c *= nu + i;

                    table.Add(c);
                }

                return table[k];
            }
        }

        private class GammaPNTable {
            private readonly MultiPrecision<N> r;
            private readonly GammaTable positive_table, negative_table;
            private readonly List<MultiPrecision<N>> table = [];

            public GammaPNTable(MultiPrecision<N> nu) {
                this.r = MultiPrecision<N>.Pow(4d, nu);
                this.positive_table = new(nu);
                this.negative_table = new(-nu);
            }

            public MultiPrecision<N> this[int k] => Value(k);

            public MultiPrecision<N> Value(int k) {
                Debug.Assert(k >= 0);

                if (k < table.Count) {
                    return table[k];
                }

                for (int i = table.Count; i <= k; i++) {
                    MultiPrecision<N> c = r * positive_table[i] / negative_table[i];

                    table.Add(c);
                }

                return table[k];
            }
        }

        private class YCoefTable {
            private MultiPrecision<N> c;
            private readonly List<MultiPrecision<N>> table = [];

            public YCoefTable() {
                this.c = 1d;
                this.table.Add(1d);
            }

            public MultiPrecision<N> this[int k] => Value(k);

            public MultiPrecision<N> Value(int k) {
                Debug.Assert(k >= 0);

                if (k < table.Count) {
                    return table[k];
                }

                for (long i = table.Count; i <= k; i++) {
                    c *= checked(32 * i * (2 * i - 1));

                    table.Add(MultiPrecision<N>.Rcp(c));
                }

                return table[k];
            }
        }

        private class Y0CoefTable {
            private readonly List<MultiPrecision<N>> table = [];

            public Y0CoefTable() {
                this.table.Add(MultiPrecision<N>.Rcp(4));
            }

            public MultiPrecision<N> this[int k] => Value(k);

            public MultiPrecision<N> Value(int k) {
                Debug.Assert(k >= 0);

                if (k < table.Count) {
                    return table[k];
                }

                for (long i = table.Count; i <= k; i++) {
                    MultiPrecision<N> c = MultiPrecision<N>.Rcp(checked(4 * (2 * i + 1) * (2 * i + 1) * (2 * i + 1)));

                    table.Add(c);
                }

                return table[k];
            }
        }

        private class YNCoefTable {
            private readonly int n;
            private readonly List<MultiPrecision<N>> table = [];

            public YNCoefTable(int n) {
                this.n = n;
            }

            public MultiPrecision<N> this[int k] => Value(k);

            public MultiPrecision<N> Value(int k) {
                Debug.Assert(k >= 0);

                if (k < table.Count) {
                    return table[k];
                }

                for (long i = table.Count; i <= k; i++) {
                    MultiPrecision<N> c = (MultiPrecision<N>)(n + 4 * i + 2) /
                        (MultiPrecision<N>)checked(4 * (2 * i + 1) * (2 * i + 1) * (n + 2 * i + 1) * (n + 2 * i + 1));

                    table.Add(c);
                }

                return table[k];
            }
        }

        private static class YNFiniteCoefTable {
            public static ReadOnlyCollection<MultiPrecision<N>> Value(int n) {
                Debug.Assert(n >= 0);

                List<MultiPrecision<N>> frac = [1], coef = [];

                for (int i = 1; i < n; i++) {
                    frac.Add(i * frac[^1]);
                }

                for (int i = 0; i < n; i++) {
                    coef.Add(frac[^(i + 1)] / frac[i]);
                }

                return new(coef);
            }
        }

        private class KCoefTable {
            private MultiPrecision<N> c;
            private readonly List<MultiPrecision<N>> table = [];

            public KCoefTable() {
                this.c = 1d;
                this.table.Add(1d);
            }

            public MultiPrecision<N> this[int k] => Value(k);

            public MultiPrecision<N> Value(int k) {
                Debug.Assert(k >= 0);

                if (k < table.Count) {
                    return table[k];
                }

                for (long i = table.Count; i <= k; i++) {
                    c *= 4 * i;

                    table.Add(MultiPrecision<N>.Rcp(c));
                }

                return table[k];
            }
        }

        private class K0CoefTable {
            private MultiPrecision<N> c;
            private readonly List<MultiPrecision<N>> table = [];

            public K0CoefTable() {
                this.c = 1d;
                this.table.Add(1d);
            }

            public MultiPrecision<N> this[int k] => Value(k);

            public MultiPrecision<N> Value(int k) {
                Debug.Assert(k >= 0);

                if (k < table.Count) {
                    return table[k];
                }

                for (long i = table.Count; i <= k; i++) {
                    c *= checked(4 * i * i);

                    table.Add(MultiPrecision<N>.Rcp(c));
                }

                return table[k];
            }
        }

        private class K1CoefTable {
            private MultiPrecision<N> c;
            private readonly List<MultiPrecision<N>> table = [];

            public K1CoefTable() {
                this.c = 1d;
                this.table.Add(1d);
            }

            public MultiPrecision<N> this[int k] => Value(k);

            public MultiPrecision<N> Value(int k) {
                Debug.Assert(k >= 0);

                if (k < table.Count) {
                    return table[k];
                }

                for (int i = table.Count; i <= k; i++) {
                    c *= checked(4 * i * (i + 1));

                    table.Add(MultiPrecision<N>.Rcp(c));
                }

                return table[k];
            }
        }
    }
}
