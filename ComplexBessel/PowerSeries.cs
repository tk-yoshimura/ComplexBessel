﻿using MultiPrecision;
using MultiPrecisionComplex;

namespace ComplexBessel {
    public class PowerSeries<N> where N : struct, IConstant {
        private static readonly Dictionary<MultiPrecision<N>, DoubleFactDenomTable> dfactdenom_coef_table = new();
        private static readonly Dictionary<MultiPrecision<N>, X2DenomTable> x2denom_coef_table = new();
        private static readonly Dictionary<MultiPrecision<N>, GammaDenomTable> gammadenom_coef_table = new();
        private static readonly Dictionary<MultiPrecision<N>, GammaTable> gamma_coef_table = new();
        private static readonly Dictionary<MultiPrecision<N>, GammaPNTable> gammapn_coef_table = new();
        private static readonly YCoefTable y_coef_table = new();
        private static readonly Y0CoefTable y0_coef_table = new();
        private static readonly Y1CoefTable y1_coef_table = new();
        private static readonly KCoefTable k_coef_table = new();
        private static readonly K0CoefTable k0_coef_table = new();
        private static readonly K1CoefTable k1_coef_table = new();

        public static Complex<N> BesselJ(MultiPrecision<N> nu, Complex<N> z) {
            if (MultiPrecision<N>.IsNegative(nu) && BesselUtil<N>.NearlyInteger(nu, out int n)) {
                Complex<N> y = BesselJ(-nu, z);

                return ((n & 1) == 0) ? y : -y;
            }
            else {
                Complex<N> y = BesselJIKernel(nu, z, sign_switch: true, terms: 256);

                return y;
            }
        }

        public static Complex<N> BesselY(MultiPrecision<N> nu, Complex<N> z) {
            if (BesselUtil<N>.NearlyInteger(nu, out int n)) {
                Complex<N> y = BesselYKernel(n, z, terms: 256);

                return y;
            }
            else if (nu < 0d && MultiPrecision<N>.Abs((nu - MultiPrecision<N>.Floor(nu)) - 0.5d) < 0.0625d) {
                Complex<N> y = BesselYKernel(nu, z, terms: 256);

                return y;
            }
            else {
                Complex<N> y = BesselYKernel(nu, z, terms: 256);

                return y;
            }
        }

        public static Complex<N> BesselI(MultiPrecision<N> nu, Complex<N> z, bool scale = false) {
            if (MultiPrecision<N>.IsNegative(nu) && BesselUtil<N>.NearlyInteger(nu, out _)) {
                Complex<N> y = BesselI(-nu, z);

                if (scale) {
                    y *= Complex<N>.Exp(-z);
                }

                return y;
            }
            else {
                Complex<N> y = BesselJIKernel(nu, z, sign_switch: false, terms: 256);

                if (scale) {
                    y *= Complex<N>.Exp(-z);
                }

                return y;
            }
        }

        public static Complex<N> BesselK(MultiPrecision<N> nu, Complex<N> z, bool scale = false) {
            if (BesselUtil<N>.NearlyInteger(nu, out int n)) {
                Complex<N> y = BesselKKernel(n, z, terms: 256);

                if (scale) {
                    y *= Complex<N>.Exp(z);
                }

                return y;
            }
            else {
                Complex<N> y = BesselKKernel(nu, z, terms: 256);

                if (scale) {
                    y *= Complex<N>.Exp(z);
                }

                return y;
            }
        }

        private static Complex<N> BesselJIKernel(MultiPrecision<N> nu, Complex<N> z, bool sign_switch, int terms) {
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

            Complex<N> z2 = z * z, z4 = z2 * z2;

            Complex<N> c = 0d, u = Complex<N>.Pow(Complex<N>.Ldexp(z, -1), nu);

            for (int k = 0, conv_times = 0; k <= terms && conv_times < 2; k++) {
                Complex<N> w = z2 * d[k];
                Complex<N> dc = u * r[k] * (sign_switch ? (1d - w) : (1d + w));

                Complex<N> c_next = c + dc;

                if (c == c_next || !Complex<N>.IsFinite(c_next)) {
                    conv_times++;
                }
                else {
                    conv_times = 0;
                }

                c = c_next;
                u *= z4;

                if (!Complex<N>.IsFinite(c)) {
                    break;
                }
            }

            return c;
        }

        private static Complex<N> BesselYKernel(MultiPrecision<N> nu, Complex<N> z, int terms) {
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

            MultiPrecision<N> cos = MultiPrecision<N>.CosPI(nu), sin = MultiPrecision<N>.SinPI(nu);
            Complex<N> p = Complex<N>.IsZero(cos) ? 0d : Complex<N>.Pow(z, Complex<N>.Ldexp(nu, 1)) * cos, s = Complex<N>.Ldexp(Complex<N>.Pow(Complex<N>.Ldexp(z, 1), nu), 2);

            Complex<N> z2 = z * z, z4 = z2 * z2;

            Complex<N> c = 0d, u = 1d / sin;

            for (int k = 0, t = 1, conv_times = 0; k <= terms && conv_times < 2; k++, t += 2) {
                Complex<N> a = t * s * g[t], q = gpn[t];
                Complex<N> pa = p / a, qa = q / a;

                Complex<N> dc = u * r[k] * ((4 * t * nu) * (pa + qa) - (z2 - (4 * t * t)) * (pa - qa));

                Complex<N> c_next = c + dc;

                if (c == c_next || !Complex<N>.IsFinite(c_next)) {
                    conv_times++;
                }
                else {
                    conv_times = 0;
                }

                c = c_next;
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

                return ((n & 1) == 0) ? y : -y;
            }
            else {
                if (n == 0) {
                    return BesselY0Kernel(z, terms);
                }
                if (n == 1) {
                    return BesselY1Kernel(z, terms);
                }
            }

            Complex<N> v = 1d / z;
            Complex<N> y0 = BesselY0Kernel(z, terms);
            Complex<N> y1 = BesselY1Kernel(z, terms);

            long exp_sum = 0;

            for (int k = 1; k < n; k++) {
                long exp = long.Max(y0.Exponent, y1.Exponent);
                (y0, y1) = (Complex<N>.Ldexp(y0, -exp), Complex<N>.Ldexp(y1, -exp));
                (y1, y0) = ((2 * k) * v * y1 - y0, y1);

                exp_sum += exp;
            }

            y1 = Complex<N>.Ldexp(y1, -exp_sum);

            return y1;
        }

        private static Complex<N> BesselY0Kernel(Complex<N> z, int terms) {
            if (!dfactdenom_coef_table.ContainsKey(0)) {
                dfactdenom_coef_table.Add(0, new DoubleFactDenomTable(0));
            }
            if (!x2denom_coef_table.ContainsKey(0)) {
                x2denom_coef_table.Add(0, new X2DenomTable(0));
            }

            DoubleFactDenomTable r = dfactdenom_coef_table[0];
            X2DenomTable d = x2denom_coef_table[0];
            Y0CoefTable q = y0_coef_table;

            Complex<N> h = Complex<N>.Log(Complex<N>.Ldexp(z, -1)) + MultiPrecision<N>.EulerGamma;

            Complex<N> z2 = z * z, z4 = z2 * z2;

            Complex<N> c = 0d, u = Complex<N>.Ldexp(MultiPrecision<N>.RcpPI, 1);

            for (int k = 0, conv_times = 0; k <= terms && conv_times < 2; k++) {
                Complex<N> dc = u * r[k] * ((h - MultiPrecision<N>.HarmonicNumber(2 * k)) * (1d - z2 * d[k]) + z2 * q[k]);

                Complex<N> c_next = c + dc;

                if (c == c_next || !Complex<N>.IsFinite(c_next)) {
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

        private static Complex<N> BesselY1Kernel(Complex<N> z, int terms) {
            if (!dfactdenom_coef_table.ContainsKey(1)) {
                dfactdenom_coef_table.Add(1, new DoubleFactDenomTable(1));
            }
            if (!x2denom_coef_table.ContainsKey(1)) {
                x2denom_coef_table.Add(1, new X2DenomTable(1));
            }

            DoubleFactDenomTable r = dfactdenom_coef_table[1];
            X2DenomTable d = x2denom_coef_table[1];
            Y1CoefTable q = y1_coef_table;

            Complex<N> h = Complex<N>.Ldexp(Complex<N>.Log(Complex<N>.Ldexp(z, -1)) + MultiPrecision<N>.EulerGamma, 1);

            Complex<N> z2 = z * z, z4 = z2 * z2;

            Complex<N> c = -2d / (z * MultiPrecision<N>.PI), u = z / Complex<N>.Ldexp(MultiPrecision<N>.PI, 1);

            for (int k = 0, conv_times = 0; k <= terms && conv_times < 2; k++) {
                Complex<N> dc = u * r[k] * ((h - MultiPrecision<N>.HarmonicNumber(2 * k) - MultiPrecision<N>.HarmonicNumber(2 * k + 1)) * (1d - z2 * d[k]) + z2 * q[k]);

                Complex<N> c_next = c + dc;

                if (c == c_next || !Complex<N>.IsFinite(c_next)) {
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

        private static Complex<N> BesselKKernel(MultiPrecision<N> nu, Complex<N> z, int terms) {
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

            Complex<N> tp = Complex<N>.Pow(Complex<N>.Ldexp(z, -1), nu), tn = 1d / tp;

            Complex<N> z2 = z * z;

            Complex<N> c = 0d, u = MultiPrecision<N>.PI / Complex<N>.Ldexp(MultiPrecision<N>.SinPI(nu), 1);

            for (int k = 0, conv_times = 0; k <= terms && conv_times < 2; k++) {
                Complex<N> dc = u * r[k] * (tn * gn[k] - tp * gp[k]);

                Complex<N> c_next = c + dc;

                if (c == c_next || !Complex<N>.IsFinite(c_next)) {
                    conv_times++;
                }
                else {
                    conv_times = 0;
                }

                c = c_next;
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
            if (n == 1) {
                return BesselK1Kernel(z, terms);
            }

            Complex<N> v = 1d / z;
            Complex<N> y0 = BesselK0Kernel(z, terms);
            Complex<N> y1 = BesselK1Kernel(z, terms);

            for (int k = 1; k < n; k++) {
                (y1, y0) = ((2 * k) * v * y1 + y0, y1);
            }

            return y1;
        }

        private static Complex<N> BesselK0Kernel(Complex<N> z, int terms) {
            K0CoefTable r = k0_coef_table;
            Complex<N> h = -Complex<N>.Log(Complex<N>.Ldexp(z, -1)) - MultiPrecision<N>.EulerGamma;

            Complex<N> z2 = z * z;

            Complex<N> c = 0d, u = 1d;

            for (int k = 0, conv_times = 0; k <= terms && conv_times < 2; k++) {
                Complex<N> dc = u * r[k] * (h + MultiPrecision<N>.HarmonicNumber(k));

                Complex<N> c_next = c + dc;

                if (c == c_next || !Complex<N>.IsFinite(c_next)) {
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

        private static Complex<N> BesselK1Kernel(Complex<N> z, int terms) {
            K1CoefTable r = k1_coef_table;
            Complex<N> h = Complex<N>.Log(Complex<N>.Ldexp(z, -1)) + MultiPrecision<N>.EulerGamma;

            Complex<N> z2 = z * z;

            Complex<N> c = 1d / z, u = Complex<N>.Ldexp(z, -1);

            for (int k = 0, conv_times = 0; k <= terms && conv_times < 2; k++) {
                Complex<N> dc = u * r[k] * (h - Complex<N>.Ldexp(MultiPrecision<N>.HarmonicNumber(k) + MultiPrecision<N>.HarmonicNumber(k + 1), -1));

                Complex<N> c_next = c + dc;

                if (c == c_next || !Complex<N>.IsFinite(c_next)) {
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
            private MultiPrecision<N> c;
            private readonly MultiPrecision<N> nu;
            private readonly List<MultiPrecision<N>> table = new();

            public DoubleFactDenomTable(MultiPrecision<N> nu) {
                this.c = MultiPrecision<N>.Gamma(nu + 1d);
                this.nu = nu;
                this.table.Add(MultiPrecision<N>.Rcp(c));
            }

            public MultiPrecision<N> this[int n] => Value(n);

            public MultiPrecision<N> Value(int n) {
                ArgumentOutOfRangeException.ThrowIfNegative(n, nameof(n));

                if (n < table.Count) {
                    return table[n];
                }

                for (int k = table.Count; k <= n; k++) {
                    c *= (nu + (2 * k)) * (nu + (2 * k - 1)) * (32 * k * (2 * k - 1));

                    table.Add(MultiPrecision<N>.Rcp(c));
                }

                return table[n];
            }
        }

        private class X2DenomTable {
            private readonly MultiPrecision<N> nu;
            private readonly List<MultiPrecision<N>> table = new();

            public X2DenomTable(MultiPrecision<N> nu) {
                MultiPrecision<N> a = MultiPrecision<N>.Rcp(4d * (nu + 1d));

                this.nu = nu;
                this.table.Add(a);
            }

            public MultiPrecision<N> this[int n] => Value(n);

            public MultiPrecision<N> Value(int n) {
                ArgumentOutOfRangeException.ThrowIfNegative(n, nameof(n));

                if (n < table.Count) {
                    return table[n];
                }

                for (int k = table.Count; k <= n; k++) {
                    MultiPrecision<N> a = MultiPrecision<N>.Rcp(4d * (2 * k + 1) * (2 * k + 1 + nu));

                    table.Add(a);
                }

                return table[n];
            }
        }

        private class GammaDenomTable {
            private MultiPrecision<N> c;
            private readonly MultiPrecision<N> nu;
            private readonly List<MultiPrecision<N>> table = new();

            public GammaDenomTable(MultiPrecision<N> nu) {
                this.c = MultiPrecision<N>.Gamma(nu + 1d);
                this.nu = nu;
                this.table.Add(MultiPrecision<N>.Rcp(c));
            }

            public MultiPrecision<N> this[int n] => Value(n);

            public MultiPrecision<N> Value(int n) {
                ArgumentOutOfRangeException.ThrowIfNegative(n, nameof(n));

                if (n < table.Count) {
                    return table[n];
                }

                for (int k = table.Count; k <= n; k++) {
                    c *= nu + k;

                    table.Add(MultiPrecision<N>.Rcp(c));
                }

                return table[n];
            }
        }

        private class GammaTable {
            private MultiPrecision<N> c;
            private readonly MultiPrecision<N> nu;
            private readonly List<MultiPrecision<N>> table = new();

            public GammaTable(MultiPrecision<N> nu) {
                this.c = MultiPrecision<N>.Gamma(nu + 1d);
                this.nu = nu;
                this.table.Add(c);
            }

            public MultiPrecision<N> this[int n] => Value(n);

            public MultiPrecision<N> Value(int n) {
                ArgumentOutOfRangeException.ThrowIfNegative(n, nameof(n));

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
            private readonly MultiPrecision<N> r;
            private readonly GammaTable positive_table, negative_table;
            private readonly List<MultiPrecision<N>> table = new();

            public GammaPNTable(MultiPrecision<N> nu) {
                this.r = MultiPrecision<N>.Pow(4, nu);
                this.positive_table = new(nu);
                this.negative_table = new(-nu);
            }

            public MultiPrecision<N> this[int n] => Value(n);

            public MultiPrecision<N> Value(int n) {
                ArgumentOutOfRangeException.ThrowIfNegative(n, nameof(n));

                if (n < table.Count) {
                    return table[n];
                }

                for (int k = table.Count; k <= n; k++) {
                    MultiPrecision<N> c = r * positive_table[k] / negative_table[k];

                    table.Add(c);
                }

                return table[n];
            }
        }

        private class YCoefTable {
            private MultiPrecision<N> c;
            private readonly List<MultiPrecision<N>> table = new();

            public YCoefTable() {
                this.c = 1d;
                this.table.Add(1d);
            }

            public MultiPrecision<N> this[int n] => Value(n);

            public MultiPrecision<N> Value(int n) {
                ArgumentOutOfRangeException.ThrowIfNegative(n, nameof(n));

                if (n < table.Count) {
                    return table[n];
                }

                for (int k = table.Count; k <= n; k++) {
                    c *= (32 * k * (2 * k - 1));

                    table.Add(MultiPrecision<N>.Rcp(c));
                }

                return table[n];
            }
        }

        private class Y0CoefTable {
            private readonly List<MultiPrecision<N>> table = new();

            public Y0CoefTable() {
                this.table.Add(MultiPrecision<N>.Rcp(4));
            }

            public MultiPrecision<N> this[int n] => Value(n);

            public MultiPrecision<N> Value(int n) {
                ArgumentOutOfRangeException.ThrowIfNegative(n, nameof(n));

                if (n < table.Count) {
                    return table[n];
                }

                for (int k = table.Count; k <= n; k++) {
                    MultiPrecision<N> c = MultiPrecision<N>.Rcp((4 * (2 * k + 1) * (2 * k + 1) * (2 * k + 1)));

                    table.Add(c);
                }

                return table[n];
            }
        }

        private class Y1CoefTable {
            private readonly List<MultiPrecision<N>> table = new();

            public Y1CoefTable() {
                this.table.Add(MultiPrecision<N>.Ldexp(3, -4));
            }

            public MultiPrecision<N> this[int n] => Value(n);

            public MultiPrecision<N> Value(int n) {
                ArgumentOutOfRangeException.ThrowIfNegative(n, nameof(n));

                if (n < table.Count) {
                    return table[n];
                }

                for (int k = table.Count; k <= n; k++) {
                    MultiPrecision<N> c = (MultiPrecision<N>)(4 * k + 3) / (MultiPrecision<N>)(4 * (2 * k + 1) * (2 * k + 1) * (2 * k + 2) * (2 * k + 2));

                    table.Add(c);
                }

                return table[n];
            }
        }

        private class KCoefTable {
            private MultiPrecision<N> c;
            private readonly List<MultiPrecision<N>> table = new();

            public KCoefTable() {
                this.c = 1d;
                this.table.Add(1d);
            }

            public MultiPrecision<N> this[int n] => Value(n);

            public MultiPrecision<N> Value(int n) {
                ArgumentOutOfRangeException.ThrowIfNegative(n, nameof(n));

                if (n < table.Count) {
                    return table[n];
                }

                for (int k = table.Count; k <= n; k++) {
                    c *= (4 * k);

                    table.Add(MultiPrecision<N>.Rcp(c));
                }

                return table[n];
            }
        }

        private class K0CoefTable {
            private MultiPrecision<N> c;
            private readonly List<MultiPrecision<N>> table = new();

            public K0CoefTable() {
                this.c = 1d;
                this.table.Add(1d);
            }

            public MultiPrecision<N> this[int n] => Value(n);

            public MultiPrecision<N> Value(int n) {
                ArgumentOutOfRangeException.ThrowIfNegative(n, nameof(n));

                if (n < table.Count) {
                    return table[n];
                }

                for (int k = table.Count; k <= n; k++) {
                    c *= (4 * k * k);

                    table.Add(MultiPrecision<N>.Rcp(c));
                }

                return table[n];
            }
        }

        private class K1CoefTable {
            private MultiPrecision<N> c;
            private readonly List<MultiPrecision<N>> table = new();

            public K1CoefTable() {
                this.c = 1d;
                this.table.Add(1d);
            }

            public MultiPrecision<N> this[int n] => Value(n);

            public MultiPrecision<N> Value(int n) {
                ArgumentOutOfRangeException.ThrowIfNegative(n, nameof(n));

                if (n < table.Count) {
                    return table[n];
                }

                for (int k = table.Count; k <= n; k++) {
                    c *= (4 * k * (k + 1));

                    table.Add(MultiPrecision<N>.Rcp(c));
                }

                return table[n];
            }
        }
    }
}