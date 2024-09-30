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

                return ((n & 1) == 0) ? y : -y;
            }
            else {
                Complex y = BesselJIKernel(nu, z, sign_switch: true, terms: 256);

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
            else if (nu < 0d && ddouble.Abs(nu - ddouble.Floor(nu) - 0.5d) < 0.0625d) {
                Complex y = BesselYKernel(nu, z, terms: 256);

                return y;
            }
            else {
                Complex y = BesselYKernel(nu, z, terms: 256);

                return y;
            }
        }

        public static Complex BesselI(ddouble nu, Complex z, bool scale = false) {
            Debug.Assert(ddouble.IsPositive(z.R));
            Debug.Assert(ddouble.IsPositive(z.I));

            if (ddouble.IsNegative(nu) && BesselUtil.NearlyInteger(nu, out _)) {
                Complex y = BesselI(-nu, z);

                if (scale) {
                    y *= Complex.Exp(-z);
                }

                return y;
            }
            else {
                Complex y = BesselJIKernel(nu, z, sign_switch: false, terms: 256);

                if (scale) {
                    y *= Complex.Exp(-z);
                }

                return y;
            }
        }

        public static Complex BesselK(ddouble nu, Complex z, bool scale = false) {
            Debug.Assert(ddouble.IsPositive(z.R));
            Debug.Assert(ddouble.IsPositive(z.I));

            if (BesselUtil.NearlyInteger(nu, out int n)) {
                Complex y = BesselKKernel(n, z, terms: 256);

                if (scale) {
                    y *= Complex.Exp(z);
                }

                return y;
            }
            else {
                Complex y = BesselKKernel(nu, z, terms: 256);

                if (scale) {
                    y *= Complex.Exp(z);
                }

                return y;
            }
        }

        private static Complex BesselJIKernel(ddouble nu, Complex z, bool sign_switch, int terms) {
            if (!dfactdenom_coef_table.TryGetValue(nu, out DoubleFactDenomTable dfactdenom_table)) {
                dfactdenom_table = new DoubleFactDenomTable(nu);
                dfactdenom_coef_table.Add(nu, dfactdenom_table);
            }
            if (!x2denom_coef_table.TryGetValue(nu, out X2DenomTable x2denom_table)) {
                x2denom_table = new X2DenomTable(nu);
                x2denom_coef_table.Add(nu, x2denom_table);
            }

            DoubleFactDenomTable r = dfactdenom_table;
            X2DenomTable d = x2denom_table;

            Complex z2 = z * z, z4 = z2 * z2;

            Complex c = 0d, u = Complex.Pow(Complex.Ldexp(z, -1), nu);

            for (int k = 0, conv_times = 0; k <= terms && conv_times < 2; k++) {
                Complex w = z2 * d[k];
                Complex dc = u * r[k] * (sign_switch ? (1d - w) : (1d + w));

                Complex c_next = c + dc;

                if (c == c_next || !Complex.IsFinite(c_next)) {
                    conv_times++;
                }
                else {
                    conv_times = 0;
                }

                c = c_next;
                u *= z4;

                if (!Complex.IsFinite(c)) {
                    break;
                }
            }

            return c;
        }

        private static Complex BesselYKernel(ddouble nu, Complex z, int terms) {
            if (!gamma_coef_table.TryGetValue(nu, out GammaTable gamma_table)) {
                gamma_table = new GammaTable(nu);
                gamma_coef_table.Add(nu, gamma_table);
            }
            if (!gammapn_coef_table.TryGetValue(nu, out GammaPNTable gammapn_table)) {
                gammapn_table = new GammaPNTable(nu);
                gammapn_coef_table.Add(nu, gammapn_table);
            }

            YCoefTable r = y_coef_table;
            GammaTable g = gamma_table;
            GammaPNTable gpn = gammapn_table;

            ddouble cos = SinCosPICache.CosPI(nu), sin = SinCosPICache.SinPI(nu);
            Complex p = Complex.IsZero(cos) ? 0d : Complex.Pow(z, Complex.Ldexp(nu, 1)) * cos, s = Complex.Ldexp(Complex.Pow(Complex.Ldexp(z, 1), nu), 2);

            Complex z2 = z * z, z4 = z2 * z2;

            Complex c = 0d, u = 1d / sin;

            for (int k = 0, t = 1, conv_times = 0; k <= terms && conv_times < 2; k++, t += 2) {
                Complex a = t * s * g[t], q = gpn[t];
                Complex pa = p / a, qa = q / a;

                Complex dc = u * r[k] * (4 * t * nu * (pa + qa) - (z2 - (4 * t * t)) * (pa - qa));

                Complex c_next = c + dc;

                if (c == c_next || !Complex.IsFinite(c_next)) {
                    conv_times++;
                }
                else {
                    conv_times = 0;
                }

                c = c_next;
                u *= z4;

                if (!Complex.IsFinite(c)) {
                    break;
                }
            }

            return c;
        }

        private static Complex BesselYKernel(int n, Complex z, int terms) {
            if (n < 0) {
                Complex y = BesselYKernel(-n, z, terms);

                return ((n & 1) == 0) ? y : -y;
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

            for (int k = 0, conv_times = 0; k <= terms && conv_times < 2; k++) {
                Complex dc = u * r[k] * ((h - ddouble.HarmonicNumber(2 * k)) * (1d - z2 * d[k]) + z2 * q[k]);

                Complex c_next = c + dc;

                if (c == c_next || !Complex.IsFinite(c_next)) {
                    conv_times++;
                }
                else {
                    conv_times = 0;
                }

                c = c_next;
                u *= z4;
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

            for (int k = 0, conv_times = 0; k <= terms && conv_times < 2; k++) {
                Complex dc = u * r[k] * ((h - ddouble.HarmonicNumber(2 * k) - ddouble.HarmonicNumber(2 * k + 1)) * (1d - z2 * d[k]) + z2 * q[k]);

                Complex c_next = c + dc;

                if (c == c_next || !Complex.IsFinite(c_next)) {
                    conv_times++;
                }
                else {
                    conv_times = 0;
                }

                c = c_next;
                u *= z4;
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

            Complex c = 0;
            Complex z2 = z * z, z4 = z2 * z2;
            Complex u = 1, v = 1, w = z2 / 4;

            for (int k = 0; k < n; k++) {
                c += v * f[k];
                v *= w;
            }
            c /= -v;

            Complex h = Complex.Ldexp(Complex.Log(Complex.Ldexp(z, -1)) + ddouble.EulerGamma, 1);

            for (int k = 0, conv_times = 0; k <= terms && conv_times < 2; k++) {
                Complex dc = u * r[k] * ((h - ddouble.HarmonicNumber(2 * k) - ddouble.HarmonicNumber(2 * k + n)) * (1d - z2 * d[k]) + z2 * q[k]);

                Complex c_next = c + dc;

                if (c == c_next || !Complex.IsFinite(c_next)) {
                    conv_times++;
                }
                else {
                    conv_times = 0;
                }

                c = c_next;
                u *= z4;
            }

            Complex y = c * ddouble.RcpPI * Complex.Pow(Complex.Ldexp(z, -1), n);

            return y;
        }

        private static Complex BesselKKernel(ddouble nu, Complex z, int terms) {
            if (!gammadenom_coef_table.TryGetValue(nu, out GammaDenomTable gammadenomp_table)) {
                gammadenomp_table = new GammaDenomTable(nu);
                gammadenom_coef_table.Add(nu, gammadenomp_table);
            }
            if (!gammadenom_coef_table.TryGetValue(-nu, out GammaDenomTable gammadenomn_table)) {
                gammadenomn_table = new GammaDenomTable(-nu);
                gammadenom_coef_table.Add(-nu, gammadenomn_table);
            }

            KCoefTable r = k_coef_table;
            GammaDenomTable gp = gammadenomp_table, gn = gammadenomn_table;

            Complex tp = Complex.Pow(Complex.Ldexp(z, -1), nu), tn = 1d / tp;

            Complex z2 = z * z;

            Complex c = 0d, u = ddouble.PI / Complex.Ldexp(ddouble.SinPI(nu), 1);

            for (int k = 0, conv_times = 0; k <= terms && conv_times < 2; k++) {
                Complex dc = u * r[k] * (tn * gn[k] - tp * gp[k]);

                Complex c_next = c + dc;

                if (c == c_next || !Complex.IsFinite(c_next)) {
                    conv_times++;
                }
                else {
                    conv_times = 0;
                }

                c = c_next;
                u *= z2;

                if (!Complex.IsFinite(c)) {
                    break;
                }
            }

            return c;
        }

        private static Complex BesselKKernel(int n, Complex z, int terms) {
            if (n == 0) {
                return BesselK0Kernel(z, terms);
            }
            if (n == 1) {
                return BesselK1Kernel(z, terms);
            }

            Complex v = 1d / z;
            Complex y0 = BesselK0Kernel(z, terms);
            Complex y1 = BesselK1Kernel(z, terms);

            for (int k = 1; k < n; k++) {
                (y1, y0) = (2 * k * v * y1 + y0, y1);
            }

            return y1;
        }

        private static Complex BesselK0Kernel(Complex z, int terms) {
            K0CoefTable r = k0_coef_table;
            Complex h = -Complex.Log(Complex.Ldexp(z, -1)) - ddouble.EulerGamma;

            Complex z2 = z * z;

            Complex c = 0d, u = 1d;

            for (int k = 0, conv_times = 0; k <= terms && conv_times < 2; k++) {
                Complex dc = u * r[k] * (h + ddouble.HarmonicNumber(k));

                Complex c_next = c + dc;

                if (c == c_next || !Complex.IsFinite(c_next)) {
                    conv_times++;
                }
                else {
                    conv_times = 0;
                }

                c = c_next;
                u *= z2;
            }

            return c;
        }

        private static Complex BesselK1Kernel(Complex z, int terms) {
            K1CoefTable r = k1_coef_table;
            Complex h = Complex.Log(Complex.Ldexp(z, -1)) + ddouble.EulerGamma;

            Complex z2 = z * z;

            Complex c = 1d / z, u = Complex.Ldexp(z, -1);

            for (int k = 0, conv_times = 0; k <= terms && conv_times < 2; k++) {
                Complex dc = u * r[k] * (h - Complex.Ldexp(ddouble.HarmonicNumber(k) + ddouble.HarmonicNumber(k + 1), -1));

                Complex c_next = c + dc;

                if (c == c_next || !Complex.IsFinite(c_next)) {
                    conv_times++;
                }
                else {
                    conv_times = 0;
                }

                c = c_next;
                u *= z2;
            }

            return c;
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

            public ddouble this[int n] => Value(n);

            public ddouble Value(int n) {
                Debug.Assert(n >= 0);

                if (n < table.Count) {
                    return table[n];
                }

                for (long k = table.Count; k <= n; k++) {
                    c *= checked((nu + (2 * k)) * (nu + (2 * k - 1)) * (32 * k * (2 * k - 1)));

                    table.Add(ddouble.Rcp(c));
                }

                return table[n];
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

            public ddouble this[int n] => Value(n);

            public ddouble Value(int n) {
                Debug.Assert(n >= 0);

                if (n < table.Count) {
                    return table[n];
                }

                for (long k = table.Count; k <= n; k++) {
                    ddouble a = ddouble.Rcp(checked(4d * (2 * k + 1) * (2 * k + 1 + nu)));

                    table.Add(a);
                }

                return table[n];
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

            public ddouble this[int n] => Value(n);

            public ddouble Value(int n) {
                Debug.Assert(n >= 0);

                if (n < table.Count) {
                    return table[n];
                }

                for (int k = table.Count; k <= n; k++) {
                    c *= nu + k;

                    table.Add(ddouble.Rcp(c));
                }

                return table[n];
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

            public ddouble this[int n] => Value(n);

            public ddouble Value(int n) {
                Debug.Assert(n >= 0);

                if (n < table.Count) {
                    return table[n];
                }

                for (int k = table.Count; k <= n; k++) {
                    c *= nu + k;

                    table.Add(c);
                }

                return table[n];
            }
        }

        private class GammaPNTable {
            private readonly ddouble r;
            private readonly GammaTable positive_table, negative_table;
            private readonly List<ddouble> table = [];

            public GammaPNTable(ddouble nu) {
                this.r = ddouble.Pow(4, nu);
                this.positive_table = new(nu);
                this.negative_table = new(-nu);
            }

            public ddouble this[int n] => Value(n);

            public ddouble Value(int n) {
                Debug.Assert(n >= 0);

                if (n < table.Count) {
                    return table[n];
                }

                for (int k = table.Count; k <= n; k++) {
                    ddouble c = r * positive_table[k] / negative_table[k];

                    table.Add(c);
                }

                return table[n];
            }
        }

        private class YCoefTable {
            private ddouble c;
            private readonly List<ddouble> table = [];

            public YCoefTable() {
                this.c = 1d;
                this.table.Add(1d);
            }

            public ddouble this[int n] => Value(n);

            public ddouble Value(int n) {
                Debug.Assert(n >= 0);

                if (n < table.Count) {
                    return table[n];
                }

                for (long k = table.Count; k <= n; k++) {
                    c *= checked(32 * k * (2 * k - 1));

                    table.Add(ddouble.Rcp(c));
                }

                return table[n];
            }
        }

        private class Y0CoefTable {
            private readonly List<ddouble> table = [];

            public Y0CoefTable() {
                this.table.Add(ddouble.Rcp(4));
            }

            public ddouble this[int n] => Value(n);

            public ddouble Value(int n) {
                Debug.Assert(n >= 0);

                if (n < table.Count) {
                    return table[n];
                }

                for (long k = table.Count; k <= n; k++) {
                    ddouble c = ddouble.Rcp(checked(4 * (2 * k + 1) * (2 * k + 1) * (2 * k + 1)));

                    table.Add(c);
                }

                return table[n];
            }
        }

        private class YNCoefTable {
            private readonly int nu;
            private readonly List<ddouble> table = [];

            public YNCoefTable(int nu) {
                this.nu = nu;
            }

            public ddouble this[int n] => Value(n);

            public ddouble Value(int n) {
                Debug.Assert(n >= 0);

                if (n < table.Count) {
                    return table[n];
                }

                for (long k = table.Count; k <= n; k++) {
                    ddouble c = (ddouble)(nu + 4 * k + 2) /
                        (ddouble)checked(4 * (2 * k + 1) * (2 * k + 1) * (nu + 2 * k + 1) * (nu + 2 * k + 1));

                    table.Add(c);
                }

                return table[n];
            }
        }

        private static class YNFiniteCoefTable {
            public static ReadOnlyCollection<ddouble> Value(int nu) {
                Debug.Assert(nu >= 0);

                List<ddouble> frac = [1], coef = [];

                for (int k = 1; k < nu; k++) {
                    frac.Add(k * frac[^1]);
                }

                for (int k = 0; k < nu; k++) {
                    coef.Add(frac[^(k + 1)] / frac[k]);
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

            public ddouble this[int n] => Value(n);

            public ddouble Value(int n) {
                Debug.Assert(n >= 0);

                if (n < table.Count) {
                    return table[n];
                }

                for (long k = table.Count; k <= n; k++) {
                    c *= 4 * k;

                    table.Add(ddouble.Rcp(c));
                }

                return table[n];
            }
        }

        private class K0CoefTable {
            private ddouble c;
            private readonly List<ddouble> table = [];

            public K0CoefTable() {
                this.c = 1d;
                this.table.Add(1d);
            }

            public ddouble this[int n] => Value(n);

            public ddouble Value(int n) {
                Debug.Assert(n >= 0);

                if (n < table.Count) {
                    return table[n];
                }

                for (long k = table.Count; k <= n; k++) {
                    c *= checked(4 * k * k);

                    table.Add(ddouble.Rcp(c));
                }

                return table[n];
            }
        }

        private class K1CoefTable {
            private ddouble c;
            private readonly List<ddouble> table = [];

            public K1CoefTable() {
                this.c = 1d;
                this.table.Add(1d);
            }

            public ddouble this[int n] => Value(n);

            public ddouble Value(int n) {
                Debug.Assert(n >= 0);

                if (n < table.Count) {
                    return table[n];
                }

                for (int k = table.Count; k <= n; k++) {
                    c *= checked(4 * k * (k + 1));

                    table.Add(ddouble.Rcp(c));
                }

                return table[n];
            }
        }
    }
}
