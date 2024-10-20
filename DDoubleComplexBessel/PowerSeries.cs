using DoubleDouble;
using DoubleDoubleComplex;
using System.Collections.ObjectModel;
using System.Diagnostics;

namespace DDoubleComplexBessel {
    public class PowerSeries {
        private static readonly Dictionary<ddouble, DoubleFactDenomTable> dfactdenom_coef_table = [];
        private static readonly Dictionary<ddouble, X2DenomTable> x2denom_coef_table = [];
        private static readonly Dictionary<ddouble, GammaDenomTable> gammadenom_coef_table = [];
        private static readonly Dictionary<ddouble, GammaTable> gamma_coef_table = [];
        private static readonly Dictionary<ddouble, GammaPNTable> gammapn_coef_table = [];
        private static readonly YCoefTable y_coef_table = new();
        private static readonly Y0CoefTable y0_coef_table = new();
        private static readonly Dictionary<int, YNCoefTable> yn_coef_table = [];
        private static readonly Dictionary<int, ReadOnlyCollection<ddouble>> yn_finitecoef_table = [];
        private static readonly KCoefTable k_coef_table = new();
        private static readonly K0CoefTable k0_coef_table = new();
        private static readonly K1CoefTable k1_coef_table = new();

        public static Complex BesselJ(ddouble nu, Complex z) {
            Debug.Assert(ddouble.IsPositive(z.R));
            Debug.Assert(ddouble.IsPositive(z.I));

            if (ddouble.IsNegative(nu) && BesselUtil.NearlyInteger(nu, out int n)) {
                Complex y = BesselJ(-nu, z);

                return (n & 1) == 0 ? y : -y;
            }
            else {
                Complex y = BesselJKernel(nu, z, terms: 256);

                return y;
            }
        }

        public static Complex BesselY(ddouble nu, Complex z) {
            Debug.Assert(ddouble.IsPositive(z.R));
            Debug.Assert(ddouble.IsPositive(z.I));

            if (BesselUtil.NearlyInteger(nu, out int n)) {
                Complex y = BesselYKernel(n, z, terms: 256);

                return y;
            }
            else {
                Complex y = BesselYKernel(nu, z, terms: 256);

                return y;
            }
        }

        public static Complex BesselI(ddouble nu, Complex z) {
            Debug.Assert(ddouble.IsPositive(z.R));
            Debug.Assert(ddouble.IsPositive(z.I));

            if (ddouble.IsNegative(nu) && BesselUtil.NearlyInteger(nu, out _)) {
                Complex y = BesselI(-nu, z);

                return y;
            }
            else {
                Complex y = BesselIKernel(nu, z, terms: 256);

                return y;
            }
        }

        public static Complex BesselK(ddouble nu, Complex z) {
            Debug.Assert(ddouble.IsPositive(z.R));
            Debug.Assert(ddouble.IsPositive(z.I));

            if (BesselUtil.NearlyInteger(nu, out int n)) {
                Complex y = BesselKKernel(n, z, terms: 256);

                return y;
            }
            else {
                Complex y = BesselKKernel(nu, z, terms: 256);

                return y;
            }
        }

        private static Complex BesselJKernel(ddouble nu, Complex z, int terms) {
            if (!dfactdenom_coef_table.TryGetValue(nu, out DoubleFactDenomTable r)) {
                r = new DoubleFactDenomTable(nu);
                dfactdenom_coef_table.Add(nu, r);
            }
            if (!x2denom_coef_table.TryGetValue(nu, out X2DenomTable d)) {
                d = new X2DenomTable(nu);
                x2denom_coef_table.Add(nu, d);
            }

            Complex z2 = z * z, z4 = z2 * z2;

            Complex c = 0d, u = Complex.Pow(Complex.Ldexp(z, -1), nu);

            for (int k = 0; k <= terms; k++) {
                c = SeriesUtil.Add(c, u * r[k], 1d, -z2 * d[k], out bool convergence);

                if (convergence) {
                    break;
                }

                u *= z4;

                if (!Complex.IsFinite(c)) {
                    break;
                }

                IterationLogger.Log("BesselJ Series", k);
            }

            return c;
        }

        private static Complex BesselYKernel(ddouble nu, Complex z, int terms) {
            if (!gamma_coef_table.TryGetValue(nu, out GammaTable g)) {
                g = new GammaTable(nu);
                gamma_coef_table.Add(nu, g);
            }
            if (!gammapn_coef_table.TryGetValue(nu, out GammaPNTable gpn)) {
                gpn = new GammaPNTable(nu);
                gammapn_coef_table.Add(nu, gpn);
            }

            YCoefTable r = y_coef_table;

            ddouble cos = SinCosPICache.CosPI(nu), sin = SinCosPICache.SinPI(nu);
            Complex p = Complex.IsZero(cos) ? 0d : Complex.Pow(z, Complex.Ldexp(nu, 1)) * cos, s = Complex.Ldexp(Complex.Pow(Complex.Ldexp(z, 1), nu), 2);

            Complex z2 = z * z, z4 = z2 * z2;

            Complex c = 0d, u = 1d / sin;

            for (int k = 0, t = 1; k <= terms; k++, t += 2) {
                Complex a = t * s * g[t], q = gpn[t];
                Complex pa = p / a, qa = q / a;

                Complex v = 4 * t * t - z2;
                c = SeriesUtil.Add(c, u * r[k], 4 * t * nu * (pa + qa), v * (pa - qa), out bool convergence);

                if (convergence && Complex.ILogB(v) >= -4) {
                    break;
                }

                u *= z4;

                if (!Complex.IsFinite(c)) {
                    break;
                }

                IterationLogger.Log("BesselY Series", k);
            }

            return c;
        }

        private static Complex BesselYKernel(int n, Complex z, int terms) {
            if (n < 0) {
                Complex y = BesselYKernel(-n, z, terms);

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

        private static Complex BesselY0Kernel(Complex z, int terms) {
            if (!dfactdenom_coef_table.TryGetValue(0, out DoubleFactDenomTable r)) {
                r = new DoubleFactDenomTable(0);
                dfactdenom_coef_table.Add(0, r);
            }
            if (!x2denom_coef_table.TryGetValue(0, out X2DenomTable d)) {
                d = new X2DenomTable(0);
                x2denom_coef_table.Add(0, d);
            }

            Y0CoefTable q = y0_coef_table;

            Complex h = Complex.Log(Complex.Ldexp(z, -1)) + ddouble.EulerGamma;

            Complex z2 = z * z, z4 = z2 * z2;

            Complex c = 0d, u = Complex.Ldexp(ddouble.RcpPI, 1);

            for (int k = 0; k <= terms; k++) {
                Complex s = u * r[k], t = h - ddouble.HarmonicNumber(2 * k);
                c = SeriesUtil.Add(c, s * t, 1d, -z2 * d[k], out bool convergence1);
                c = SeriesUtil.Add(c, s, z2 * q[k], out bool convergence2);

                if (convergence1 && convergence2 && Complex.ILogB(t) >= -4) {
                    break;
                }

                u *= z4;

                IterationLogger.Log("BesselY0 Series", k);
            }

            return c;
        }

        private static Complex BesselY1Kernel(Complex z, int terms) {
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

            Complex h = Complex.Ldexp(Complex.Log(Complex.Ldexp(z, -1)) + ddouble.EulerGamma, 1);

            Complex z2 = z * z, z4 = z2 * z2;

            Complex c = -2d / (z * ddouble.PI), u = z / Complex.Ldexp(ddouble.PI, 1);

            for (int k = 0; k <= terms; k++) {
                Complex s = u * r[k], t = h - ddouble.HarmonicNumber(2 * k) - ddouble.HarmonicNumber(2 * k + 1);
                c = SeriesUtil.Add(c, s * t, 1d, -z2 * d[k], out bool convergence1);
                c = SeriesUtil.Add(c, s, z2 * q[k], out bool convergence2);

                if (convergence1 && convergence2 && Complex.ILogB(t) >= -4) {
                    break;
                }

                u *= z4;

                IterationLogger.Log("BesselY1 Series", k);
            }

            return c;
        }

        private static Complex BesselYNKernel(int n, Complex z, int terms) {
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
            if (!yn_finitecoef_table.TryGetValue(n, out ReadOnlyCollection<ddouble> f)) {
                f = YNFiniteCoefTable.Value(n);
                yn_finitecoef_table.Add(n, f);
            }

            Complex c = 0d;
            Complex z2 = z * z, z4 = z2 * z2;
            Complex u = 1d, v = 1d, w = Complex.Ldexp(z2, -2);

            for (int k = 0; k < n; k++) {
                c += v * f[k];
                v *= w;
            }
            c /= -v;

            Complex h = Complex.Ldexp(Complex.Log(Complex.Ldexp(z, -1)) + ddouble.EulerGamma, 1);

            for (int k = 0; k <= terms; k++) {
                Complex s = u * r[k], t = (h - ddouble.HarmonicNumber(2 * k) - ddouble.HarmonicNumber(2 * k + n));
                c = SeriesUtil.Add(c, s * t, 1d, -z2 * d[k], out bool convergence1);
                c = SeriesUtil.Add(c, s, z2 * q[k], out bool convergence2);

                if (convergence1 && convergence2 && Complex.ILogB(t) >= -4) {
                    break;
                }

                u *= z4;

                IterationLogger.Log("BesselYN Series", k);
            }

            Complex y = c * ddouble.RcpPI * Complex.Pow(Complex.Ldexp(z, -1), n);

            return y;
        }

        private static Complex BesselIKernel(ddouble nu, Complex z, int terms) {
            if (!dfactdenom_coef_table.TryGetValue(nu, out DoubleFactDenomTable r)) {
                r = new DoubleFactDenomTable(nu);
                dfactdenom_coef_table.Add(nu, r);
            }
            if (!x2denom_coef_table.TryGetValue(nu, out X2DenomTable d)) {
                d = new X2DenomTable(nu);
                x2denom_coef_table.Add(nu, d);
            }

            Complex z2 = z * z, z4 = z2 * z2;

            Complex c = 0d, u = Complex.Pow(Complex.Ldexp(z, -1), nu);

            for (int k = 0; k <= terms; k++) {
                c = SeriesUtil.Add(c, u * r[k], 1d, z2 * d[k], out bool convergence);

                if (convergence) {
                    break;
                }

                u *= z4;

                if (!Complex.IsFinite(c)) {
                    break;
                }

                IterationLogger.Log("BesselI Series", k);
            }

            return c;
        }

        private static Complex BesselKKernel(ddouble nu, Complex z, int terms) {
            if (!gammadenom_coef_table.TryGetValue(nu, out GammaDenomTable gp)) {
                gp = new GammaDenomTable(nu);
                gammadenom_coef_table.Add(nu, gp);
            }
            if (!gammadenom_coef_table.TryGetValue(-nu, out GammaDenomTable gn)) {
                gn = new GammaDenomTable(-nu);
                gammadenom_coef_table.Add(-nu, gn);
            }

            KCoefTable r = k_coef_table;

            Complex tp = Complex.Pow(Complex.Ldexp(z, -1), nu), tn = 1d / tp;

            Complex z2 = z * z;

            Complex c = 0d, u = ddouble.PI / Complex.Ldexp(ddouble.SinPI(nu), 1);

            for (int k = 0; k <= terms; k++) {
                c = SeriesUtil.Add(c, u * r[k], tn * gn[k], -tp * gp[k], out bool convergence);

                if (convergence) {
                    break;
                }

                u *= z2;

                if (!Complex.IsFinite(c)) {
                    break;
                }

                IterationLogger.Log("BesselK Series", k);
            }

            return c;
        }

        private static Complex BesselKKernel(int n, Complex z, int terms) {
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

        private static Complex BesselK0Kernel(Complex z, int terms) {
            K0CoefTable r = k0_coef_table;
            Complex h = -Complex.Log(Complex.Ldexp(z, -1)) - ddouble.EulerGamma;

            Complex z2 = z * z;

            Complex c = 0d, u = 1d;

            for (int k = 0; k <= terms; k++) {
                c = SeriesUtil.Add(c, u * r[k], h, ddouble.HarmonicNumber(k), out bool convergence);

                if (convergence) {
                    break;
                }

                u *= z2;

                IterationLogger.Log("BesselK0 Series", k);
            }

            return c;
        }

        private static Complex BesselK1Kernel(Complex z, int terms) {
            K1CoefTable r = k1_coef_table;
            Complex h = Complex.Log(Complex.Ldexp(z, -1)) + ddouble.EulerGamma;

            Complex z2 = z * z;

            Complex c = 1d / z, u = Complex.Ldexp(z, -1);

            for (int k = 0; k <= terms; k++) {
                c = SeriesUtil.Add(c, u * r[k], h, -ddouble.Ldexp(ddouble.HarmonicNumber(k) + ddouble.HarmonicNumber(k + 1), -1), out bool convergence);

                if (convergence) {
                    break;
                }

                u *= z2;

                IterationLogger.Log("BesselK1 Series", k);
            }

            return c;
        }

        private static Complex BesselKNKernel(int n, Complex z, int terms) {
            Complex v = 1d / z;
            Complex y0 = BesselK0Kernel(z, terms);
            Complex y1 = BesselK1Kernel(z, terms);

            for (int k = 1; k < n; k++) {
                (y1, y0) = (2 * k * v * y1 + y0, y1);
            }

            return y1;
        }

        private class DoubleFactDenomTable {
            private ddouble c;
            private readonly ddouble nu;
            private readonly List<ddouble> table = [];

            public DoubleFactDenomTable(ddouble nu) {
                this.c = ddouble.Gamma(nu + 1d);
                this.nu = nu;
                this.table.Add(ddouble.Rcp(c));
            }

            public ddouble this[int k] => Value(k);

            public ddouble Value(int k) {
                Debug.Assert(k >= 0);

                if (k < table.Count) {
                    return table[k];
                }

                for (long i = table.Count; i <= k; i++) {
                    c *= checked((nu + 2 * i) * (nu + (2 * i - 1)) * (32 * i * (2 * i - 1)));

                    table.Add(ddouble.Rcp(c));
                }

                return table[k];
            }
        }

        private class X2DenomTable {
            private readonly ddouble nu;
            private readonly List<ddouble> table = [];

            public X2DenomTable(ddouble nu) {
                ddouble a = ddouble.Rcp(4d * (nu + 1d));

                this.nu = nu;
                this.table.Add(a);
            }

            public ddouble this[int k] => Value(k);

            public ddouble Value(int k) {
                Debug.Assert(k >= 0);

                if (k < table.Count) {
                    return table[k];
                }

                for (long i = table.Count; i <= k; i++) {
                    ddouble a = ddouble.Rcp(checked(4d * (2 * i + 1) * (2 * i + 1 + nu)));

                    table.Add(a);
                }

                return table[k];
            }
        }

        private class GammaDenomTable {
            private ddouble c;
            private readonly ddouble nu;
            private readonly List<ddouble> table = [];

            public GammaDenomTable(ddouble nu) {
                this.c = ddouble.Gamma(nu + 1d);
                this.nu = nu;
                this.table.Add(ddouble.Rcp(c));
            }

            public ddouble this[int k] => Value(k);

            public ddouble Value(int k) {
                Debug.Assert(k >= 0);

                if (k < table.Count) {
                    return table[k];
                }

                for (int i = table.Count; i <= k; i++) {
                    c *= nu + i;

                    table.Add(ddouble.Rcp(c));
                }

                return table[k];
            }
        }

        private class GammaTable {
            private ddouble c;
            private readonly ddouble nu;
            private readonly List<ddouble> table = [];

            public GammaTable(ddouble nu) {
                this.c = ddouble.Gamma(nu + 1d);
                this.nu = nu;
                this.table.Add(c);
            }

            public ddouble this[int k] => Value(k);

            public ddouble Value(int k) {
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
            private readonly ddouble r;
            private readonly GammaTable positive_table, negative_table;
            private readonly List<ddouble> table = [];

            public GammaPNTable(ddouble nu) {
                this.r = ddouble.Pow(4d, nu);
                this.positive_table = new(nu);
                this.negative_table = new(-nu);
            }

            public ddouble this[int k] => Value(k);

            public ddouble Value(int k) {
                Debug.Assert(k >= 0);

                if (k < table.Count) {
                    return table[k];
                }

                for (int i = table.Count; i <= k; i++) {
                    ddouble c = r * positive_table[i] / negative_table[i];

                    table.Add(c);
                }

                return table[k];
            }
        }

        private class YCoefTable {
            private ddouble c;
            private readonly List<ddouble> table = [];

            public YCoefTable() {
                this.c = 1d;
                this.table.Add(1d);
            }

            public ddouble this[int k] => Value(k);

            public ddouble Value(int k) {
                Debug.Assert(k >= 0);

                if (k < table.Count) {
                    return table[k];
                }

                for (long i = table.Count; i <= k; i++) {
                    c *= checked(32 * i * (2 * i - 1));

                    table.Add(ddouble.Rcp(c));
                }

                return table[k];
            }
        }

        private class Y0CoefTable {
            private readonly List<ddouble> table = [];

            public Y0CoefTable() {
                this.table.Add(ddouble.Rcp(4));
            }

            public ddouble this[int k] => Value(k);

            public ddouble Value(int k) {
                Debug.Assert(k >= 0);

                if (k < table.Count) {
                    return table[k];
                }

                for (long i = table.Count; i <= k; i++) {
                    ddouble c = ddouble.Rcp(checked(4 * (2 * i + 1) * (2 * i + 1) * (2 * i + 1)));

                    table.Add(c);
                }

                return table[k];
            }
        }

        private class YNCoefTable {
            private readonly int n;
            private readonly List<ddouble> table = [];

            public YNCoefTable(int n) {
                this.n = n;
            }

            public ddouble this[int k] => Value(k);

            public ddouble Value(int k) {
                Debug.Assert(k >= 0);

                if (k < table.Count) {
                    return table[k];
                }

                for (long i = table.Count; i <= k; i++) {
                    ddouble c = (ddouble)(n + 4 * i + 2) /
                        (ddouble)checked(4 * (2 * i + 1) * (2 * i + 1) * (n + 2 * i + 1) * (n + 2 * i + 1));

                    table.Add(c);
                }

                return table[k];
            }
        }

        private static class YNFiniteCoefTable {
            public static ReadOnlyCollection<ddouble> Value(int n) {
                Debug.Assert(n >= 0);

                List<ddouble> frac = [1], coef = [];

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
            private ddouble c;
            private readonly List<ddouble> table = [];

            public KCoefTable() {
                this.c = 1d;
                this.table.Add(1d);
            }

            public ddouble this[int k] => Value(k);

            public ddouble Value(int k) {
                Debug.Assert(k >= 0);

                if (k < table.Count) {
                    return table[k];
                }

                for (long i = table.Count; i <= k; i++) {
                    c *= 4 * i;

                    table.Add(ddouble.Rcp(c));
                }

                return table[k];
            }
        }

        private class K0CoefTable {
            private ddouble c;
            private readonly List<ddouble> table = [];

            public K0CoefTable() {
                this.c = 1d;
                this.table.Add(1d);
            }

            public ddouble this[int k] => Value(k);

            public ddouble Value(int k) {
                Debug.Assert(k >= 0);

                if (k < table.Count) {
                    return table[k];
                }

                for (long i = table.Count; i <= k; i++) {
                    c *= checked(4 * i * i);

                    table.Add(ddouble.Rcp(c));
                }

                return table[k];
            }
        }

        private class K1CoefTable {
            private ddouble c;
            private readonly List<ddouble> table = [];

            public K1CoefTable() {
                this.c = 1d;
                this.table.Add(1d);
            }

            public ddouble this[int k] => Value(k);

            public ddouble Value(int k) {
                Debug.Assert(k >= 0);

                if (k < table.Count) {
                    return table[k];
                }

                for (int i = table.Count; i <= k; i++) {
                    c *= checked(4 * i * (i + 1));

                    table.Add(ddouble.Rcp(c));
                }

                return table[k];
            }
        }
    }
}
