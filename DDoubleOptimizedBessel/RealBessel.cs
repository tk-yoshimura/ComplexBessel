﻿using DoubleDouble;
using System.Collections.ObjectModel;
using System.Diagnostics;
using static DDoubleOptimizedBessel.RealBesselUtil;

namespace DDoubleOptimizedBessel {

    public static class RealBessel {

        public static ddouble BesselJ(ddouble nu, ddouble x) {
            CheckNu(nu);

            if (ddouble.IsNegative(x) || ddouble.IsNaN(x)) {
                return ddouble.NaN;
            }

            if (x >= HankelThreshold) {
                return Limit.BesselJ(nu, x);
            }
            else if (x <= PowerSeriesThreshold(nu)) {
                return PowerSeries.BesselJ(nu, x);
            }
            else {
                return MillerBackward.BesselJ(nu, x);
            }
        }

        public static ddouble BesselY(ddouble nu, ddouble x) {
            CheckNu(nu);

            if (ddouble.IsNegative(x) || ddouble.IsNaN(x)) {
                return ddouble.NaN;
            }

            if (x >= HankelThreshold) {
                return Limit.BesselY(nu, x);
            }
            else if (x <= PowerSeriesThreshold(nu) - BesselJYPowerseriesBias) {
                if (NearlyInteger(nu, out _) || ddouble.Abs(ddouble.Round(nu) - nu) >= InterpolationThreshold) {
                    return PowerSeries.BesselY(nu, x);
                }
                else {
                    return CubicInterpolate.BesselYPowerSeries(nu, x);
                }
            }
            else {
                return MillerBackward.BesselY(nu, x);
            }
        }

        public static ddouble BesselI(ddouble nu, ddouble x, bool scale = false) {
            CheckNu(nu);

            if (ddouble.IsNegative(x) || ddouble.IsNaN(x)) {
                return ddouble.NaN;
            }

            if (x >= HankelThreshold) {
                return Limit.BesselI(nu, x, scale);
            }
            else {
                return PowerSeries.BesselI(nu, x, scale);
            }
        }

        public static ddouble BesselK(ddouble nu, ddouble x, bool scale = false) {
            CheckNu(nu);

            if (ddouble.IsNegative(x) || ddouble.IsNaN(x)) {
                return ddouble.NaN;
            }

            nu = ddouble.Abs(nu);

            if (x >= HankelThreshold) {
                return Limit.BesselK(nu, x, scale);
            }
            else if (x <= BesselKNearZeroThreshold) {
                if (NearlyInteger(nu, out _) || ddouble.Abs(ddouble.Round(nu) - nu) >= InterpolationThreshold) {
                    return PowerSeries.BesselK(nu, x, scale);
                }
                else {
                    return CubicInterpolate.BesselKPowerSeries(nu, x, scale);
                }
            }
            else {
                return YoshidaPade.BesselK(nu, x, scale);
            }
        }
    }

    static class RealBesselUtil {
        public const int MaxN = 16;
        public static readonly double Eps = double.ScaleB(1, -105);
        public static readonly double ExtremelyNearZero = double.ScaleB(1, -28);
        public static readonly double InterpolationThreshold = double.ScaleB(1, -25);
        public const double HankelThreshold = 38.875, MillerBackwardThreshold = 6;
        public const double BesselKPadeThreshold = 1, BesselKNearZeroThreshold = 2, BesselJYPowerseriesBias = 2;

        public static ddouble PowerSeriesThreshold(ddouble nu) {
            ddouble nu_abs = ddouble.Abs(nu);

            return 7.5 + nu_abs * (3.57e-1 + nu_abs * 5.23e-3);
        }

        public static void CheckNu(ddouble nu) {
            if (!(ddouble.Abs(nu) <= MaxN)) {
                throw new ArgumentOutOfRangeException(
                    nameof(nu),
                    $"In the calculation of the Bessel function, nu with an absolute value greater than {MaxN} is not supported."
                );
            }
        }

        public static void CheckN(int n) {
            if (n < -MaxN || n > MaxN) {
                throw new ArgumentOutOfRangeException(
                    nameof(n),
                    $"In the calculation of the Bessel function, n with an absolute value greater than {MaxN} is not supported."
                );
            }
        }

        public static bool NearlyInteger(ddouble nu, out int n) {
            n = (int)ddouble.Round(nu);

            return ddouble.Abs(nu - n) < Eps;
        }

        public static class SinCosPICache {
            private static readonly Dictionary<ddouble, ddouble> cospi_table = [];
            private static readonly Dictionary<ddouble, ddouble> sinpi_table = [];

            public static ddouble CosPI(ddouble theta) {
                if (!cospi_table.TryGetValue(theta, out ddouble cospi)) {
                    cospi = ddouble.CosPI(theta);
                    cospi_table[theta] = cospi;
                }

                return cospi;
            }

            public static ddouble SinPI(ddouble theta) {
                if (!sinpi_table.TryGetValue(theta, out ddouble sinpi)) {
                    sinpi = ddouble.SinPI(theta);
                    sinpi_table[theta] = sinpi;
                }

                return sinpi;
            }
        }

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

            public static ddouble BesselJ(ddouble nu, ddouble x) {
                Debug.Assert(ddouble.IsPositive(x));

                if (ddouble.IsNegative(nu) && NearlyInteger(nu, out int n)) {
                    ddouble y = BesselJ(-nu, x);

                    return (n & 1) == 0 ? y : -y;
                }
                else {
                    ddouble y = BesselJKernel(nu, x, terms: 27);

                    return y;
                }
            }

            public static ddouble BesselY(ddouble nu, ddouble x) {
                Debug.Assert(ddouble.IsPositive(x));

                if (NearlyInteger(nu, out int n)) {
                    ddouble y = BesselYKernel(n, x, terms: 18);

                    return y;
                }
                else {
                    ddouble y = BesselYKernel(nu, x, terms: 25);

                    return y;
                }
            }

            public static ddouble BesselI(ddouble nu, ddouble x, bool scale = false) {
                Debug.Assert(ddouble.IsPositive(x));

                if (ddouble.IsNegative(nu) && NearlyInteger(nu, out _)) {
                    ddouble y = BesselI(-nu, x);

                    if (scale) {
                        y *= ddouble.Exp(-x);
                    }

                    return y;
                }
                else {
                    ddouble y = BesselIKernel(nu, x, terms: 40);

                    if (scale) {
                        y *= ddouble.Exp(-x);
                    }

                    return y;
                }
            }

            public static ddouble BesselK(ddouble nu, ddouble x, bool scale = false) {
                Debug.Assert(ddouble.IsPositive(x));

                if (NearlyInteger(nu, out int n)) {
                    ddouble y = BesselKKernel(n, x, terms: 27);

                    if (scale) {
                        y *= ddouble.Exp(x);
                    }

                    return y;
                }
                else {
                    ddouble y = BesselKKernel(nu, x, terms: 30);

                    if (scale) {
                        y *= ddouble.Exp(x);
                    }

                    return y;
                }
            }

            private static ddouble BesselJKernel(ddouble nu, ddouble x, int terms) {
                if (!dfactdenom_coef_table.TryGetValue(nu, out DoubleFactDenomTable r)) {
                    r = new DoubleFactDenomTable(nu);
                    dfactdenom_coef_table.Add(nu, r);
                }
                if (!x2denom_coef_table.TryGetValue(nu, out X2DenomTable d)) {
                    d = new X2DenomTable(nu);
                    x2denom_coef_table.Add(nu, d);
                }

                ddouble x2 = x * x, x4 = x2 * x2;

                ddouble c = 0d, u = ddouble.Pow(ddouble.Ldexp(x, -1), nu);

                for (int k = 0, conv_times = 0; k <= terms && conv_times < 2; k++) {
                    ddouble w = x2 * d[k];
                    ddouble dc = u * r[k] * (1d - w);

                    ddouble c_next = c + dc;

                    if (c == c_next || !ddouble.IsFinite(c_next)) {
                        conv_times++;
                    }
                    else {
                        conv_times = 0;
                    }

                    c = c_next;
                    u *= x4;

                    if (!ddouble.IsFinite(c)) {
                        break;
                    }
                }

                return c;
            }

            private static ddouble BesselYKernel(ddouble nu, ddouble x, int terms) {
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
                ddouble p = ddouble.IsZero(cos) ? 0d : ddouble.Pow(x, ddouble.Ldexp(nu, 1)) * cos, s = ddouble.Ldexp(ddouble.Pow(ddouble.Ldexp(x, 1), nu), 2);

                ddouble x2 = x * x, x4 = x2 * x2;

                ddouble c = 0d, u = 1d / sin;

                for (int k = 0, t = 1, conv_times = 0; k <= terms && conv_times < 2; k++, t += 2) {
                    ddouble a = t * s * g[t], q = gpn[t];
                    ddouble pa = p / a, qa = q / a;

                    ddouble dc = u * r[k] * (4 * t * nu * (pa + qa) - (x2 - 4 * t * t) * (pa - qa));

                    ddouble c_next = c + dc;

                    if (c == c_next || !ddouble.IsFinite(c_next)) {
                        conv_times++;
                    }
                    else {
                        conv_times = 0;
                    }

                    c = c_next;
                    u *= x4;

                    if (!ddouble.IsFinite(c)) {
                        break;
                    }
                }

                return c;
            }

            private static ddouble BesselYKernel(int n, ddouble x, int terms) {
                if (n < 0) {
                    ddouble y = BesselYKernel(-n, x, terms);

                    return (n & 1) == 0 ? y : -y;
                }
                else if (n == 0) {
                    return BesselY0Kernel(x, terms);
                }
                else if (n == 1) {
                    return BesselY1Kernel(x, terms);
                }
                else {
                    return BesselYNKernel(n, x, terms);
                }
            }

            private static ddouble BesselY0Kernel(ddouble x, int terms) {
                if (!dfactdenom_coef_table.TryGetValue(0, out DoubleFactDenomTable r)) {
                    r = new DoubleFactDenomTable(0);
                    dfactdenom_coef_table.Add(0, r);
                }
                if (!x2denom_coef_table.TryGetValue(0, out X2DenomTable d)) {
                    d = new X2DenomTable(0);
                    x2denom_coef_table.Add(0, d);
                }

                Y0CoefTable q = y0_coef_table;

                ddouble h = ddouble.Log(ddouble.Ldexp(x, -1)) + ddouble.EulerGamma;

                ddouble x2 = x * x, x4 = x2 * x2;

                ddouble c = 0d, u = ddouble.Ldexp(ddouble.RcpPI, 1);

                for (int k = 0, conv_times = 0; k <= terms && conv_times < 2; k++) {
                    ddouble dc = u * r[k] * ((h - ddouble.HarmonicNumber(2 * k)) * (1d - x2 * d[k]) + x2 * q[k]);

                    ddouble c_next = c + dc;

                    if (c == c_next || !ddouble.IsFinite(c_next)) {
                        conv_times++;
                    }
                    else {
                        conv_times = 0;
                    }

                    c = c_next;
                    u *= x4;
                }

                return c;
            }

            private static ddouble BesselY1Kernel(ddouble x, int terms) {
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

                ddouble h = ddouble.Ldexp(ddouble.Log(ddouble.Ldexp(x, -1)) + ddouble.EulerGamma, 1);

                ddouble x2 = x * x, x4 = x2 * x2;

                ddouble c = -2d / (x * ddouble.PI), u = x / ddouble.Ldexp(ddouble.PI, 1);

                for (int k = 0, conv_times = 0; k <= terms && conv_times < 2; k++) {
                    ddouble dc = u * r[k] * ((h - ddouble.HarmonicNumber(2 * k) - ddouble.HarmonicNumber(2 * k + 1)) * (1d - x2 * d[k]) + x2 * q[k]);

                    ddouble c_next = c + dc;

                    if (c == c_next || !ddouble.IsFinite(c_next)) {
                        conv_times++;
                    }
                    else {
                        conv_times = 0;
                    }

                    c = c_next;
                    u *= x4;
                }

                return c;
            }

            private static ddouble BesselYNKernel(int n, ddouble x, int terms) {
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

                ddouble c = 0d;
                ddouble x2 = x * x, x4 = x2 * x2;
                ddouble u = 1d, v = 1d, w = ddouble.Ldexp(x2, -2);

                for (int k = 0; k < n; k++) {
                    c += v * f[k];
                    v *= w;
                }
                c /= -v;

                ddouble h = ddouble.Ldexp(ddouble.Log(ddouble.Ldexp(x, -1)) + ddouble.EulerGamma, 1);

                for (int k = 0, conv_times = 0; k <= terms && conv_times < 2; k++) {
                    ddouble dc = u * r[k] * ((h - ddouble.HarmonicNumber(2 * k) - ddouble.HarmonicNumber(2 * k + n)) * (1d - x2 * d[k]) + x2 * q[k]);

                    ddouble c_next = c + dc;

                    if (c == c_next || !ddouble.IsFinite(c_next)) {
                        conv_times++;
                    }
                    else {
                        conv_times = 0;
                    }

                    c = c_next;
                    u *= x4;
                }

                ddouble y = c * ddouble.RcpPI * ddouble.Pow(ddouble.Ldexp(x, -1), n);

                return y;
            }

            private static ddouble BesselIKernel(ddouble nu, ddouble x, int terms) {
                if (!dfactdenom_coef_table.TryGetValue(nu, out DoubleFactDenomTable r)) {
                    r = new DoubleFactDenomTable(nu);
                    dfactdenom_coef_table.Add(nu, r);
                }
                if (!x2denom_coef_table.TryGetValue(nu, out X2DenomTable d)) {
                    d = new X2DenomTable(nu);
                    x2denom_coef_table.Add(nu, d);
                }

                ddouble x2 = x * x, x4 = x2 * x2;

                ddouble c = 0d, u = ddouble.Pow(ddouble.Ldexp(x, -1), nu);

                for (int k = 0, conv_times = 0; k <= terms && conv_times < 2; k++) {
                    ddouble w = x2 * d[k];
                    ddouble dc = u * r[k] * (1d + w);

                    ddouble c_next = c + dc;

                    if (c == c_next || !ddouble.IsFinite(c_next)) {
                        conv_times++;
                    }
                    else {
                        conv_times = 0;
                    }

                    c = c_next;
                    u *= x4;

                    if (!ddouble.IsFinite(c)) {
                        break;
                    }
                }

                return c;
            }

            private static ddouble BesselKKernel(ddouble nu, ddouble x, int terms) {
                if (!gammadenom_coef_table.TryGetValue(nu, out GammaDenomTable gp)) {
                    gp = new GammaDenomTable(nu);
                    gammadenom_coef_table.Add(nu, gp);
                }
                if (!gammadenom_coef_table.TryGetValue(-nu, out GammaDenomTable gn)) {
                    gn = new GammaDenomTable(-nu);
                    gammadenom_coef_table.Add(-nu, gn);
                }

                KCoefTable r = k_coef_table;

                ddouble tp = ddouble.Pow(ddouble.Ldexp(x, -1), nu), tn = 1d / tp;

                ddouble x2 = x * x;

                ddouble c = 0d, u = ddouble.PI / ddouble.Ldexp(ddouble.SinPI(nu), 1);

                for (int k = 0, conv_times = 0; k <= terms && conv_times < 2; k++) {
                    ddouble dc = u * r[k] * (tn * gn[k] - tp * gp[k]);

                    ddouble c_next = c + dc;

                    if (c == c_next || !ddouble.IsFinite(c_next)) {
                        conv_times++;
                    }
                    else {
                        conv_times = 0;
                    }

                    c = c_next;
                    u *= x2;

                    if (!ddouble.IsFinite(c)) {
                        break;
                    }
                }

                return c;
            }

            private static ddouble BesselKKernel(int n, ddouble x, int terms) {
                if (n == 0) {
                    return BesselK0Kernel(x, terms);
                }
                else if (n == 1) {
                    return BesselK1Kernel(x, terms);
                }
                else {
                    return BesselKNKernel(n, x, terms);
                }
            }

            private static ddouble BesselK0Kernel(ddouble x, int terms) {
                K0CoefTable r = k0_coef_table;
                ddouble h = -ddouble.Log(ddouble.Ldexp(x, -1)) - ddouble.EulerGamma;

                ddouble x2 = x * x;

                ddouble c = 0d, u = 1d;

                for (int k = 0, conv_times = 0; k <= terms && conv_times < 2; k++) {
                    ddouble dc = u * r[k] * (h + ddouble.HarmonicNumber(k));

                    ddouble c_next = c + dc;

                    if (c == c_next || !ddouble.IsFinite(c_next)) {
                        conv_times++;
                    }
                    else {
                        conv_times = 0;
                    }

                    c = c_next;
                    u *= x2;
                }

                return c;
            }

            private static ddouble BesselK1Kernel(ddouble x, int terms) {
                K1CoefTable r = k1_coef_table;
                ddouble h = ddouble.Log(ddouble.Ldexp(x, -1)) + ddouble.EulerGamma;

                ddouble x2 = x * x;

                ddouble c = 1d / x, u = ddouble.Ldexp(x, -1);

                for (int k = 0, conv_times = 0; k <= terms && conv_times < 2; k++) {
                    ddouble dc = u * r[k] * (h - ddouble.Ldexp(ddouble.HarmonicNumber(k) + ddouble.HarmonicNumber(k + 1), -1));

                    ddouble c_next = c + dc;

                    if (c == c_next || !ddouble.IsFinite(c_next)) {
                        conv_times++;
                    }
                    else {
                        conv_times = 0;
                    }

                    c = c_next;
                    u *= x2;
                }

                return c;
            }

            private static ddouble BesselKNKernel(int n, ddouble x, int terms) {
                ddouble v = 1d / x;
                ddouble y0 = BesselK0Kernel(x, terms);
                ddouble y1 = BesselK1Kernel(x, terms);

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

        public static class Limit {
            static readonly Dictionary<ddouble, HankelExpansion> table = [];

            public static ddouble BesselJ(ddouble nu, ddouble x) {
                Debug.Assert(ddouble.IsPositive(x));

                if (!table.TryGetValue(nu, out HankelExpansion hankel)) {
                    hankel = new HankelExpansion(nu);
                    table.Add(nu, hankel);
                }

                (ddouble c_even, ddouble c_odd) = hankel.BesselJYCoef(x);

                ddouble omega = hankel.Omega(x);

                ddouble cos = ddouble.Cos(omega), sin = ddouble.Sin(omega);

                ddouble y = ddouble.Sqrt(2d / (ddouble.PI * x)) * (cos * c_even - sin * c_odd);

                return y;
            }

            public static ddouble BesselY(ddouble nu, ddouble x) {
                Debug.Assert(ddouble.IsPositive(x));

                if (!table.TryGetValue(nu, out HankelExpansion hankel)) {
                    hankel = new HankelExpansion(nu);
                    table.Add(nu, hankel);
                }

                (ddouble c_even, ddouble c_odd) = hankel.BesselJYCoef(x);

                ddouble omega = hankel.Omega(x);

                ddouble cos = ddouble.Cos(omega), sin = ddouble.Sin(omega);

                ddouble y = ddouble.Sqrt(2d / (ddouble.PI * x)) * (sin * c_even + cos * c_odd);

                return y;
            }

            public static ddouble BesselI(ddouble nu, ddouble x, bool scale = false) {
                Debug.Assert(ddouble.IsPositive(x));

                if (!table.TryGetValue(nu, out HankelExpansion hankel)) {
                    hankel = new HankelExpansion(nu);
                    table.Add(nu, hankel);
                }

                ddouble c = hankel.BesselICoef(x);

                ddouble y = ddouble.Sqrt(1d / (2d * ddouble.PI * x)) * c;

                if (!scale) {
                    y *= ddouble.Exp(x);
                }

                return y;
            }

            public static ddouble BesselK(ddouble nu, ddouble x, bool scale = false) {
                Debug.Assert(ddouble.IsPositive(nu));
                Debug.Assert(ddouble.IsPositive(x));

                if (!table.TryGetValue(nu, out HankelExpansion hankel)) {
                    hankel = new HankelExpansion(nu);
                    table.Add(nu, hankel);
                }

                ddouble c = hankel.BesselKCoef(x);

                ddouble y = ddouble.Sqrt(ddouble.PI / (2d * x)) * c;

                if (!scale) {
                    y *= ddouble.Exp(-x);
                }

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
                    for (int k = a_coef.Count; k <= n; k++) {
                        ddouble a = a_coef.Last() * (4d * Nu * Nu - checked((2 * k - 1) * (2 * k - 1))) / (k * 8);
                        a_coef.Add(a);
                    }

                    return a_coef[n];
                }

                public ddouble Omega(ddouble x) {
                    ddouble omega = x - ddouble.Ldexp(2 * Nu + 1, -2) * ddouble.PI;

                    return omega;
                }

                public (ddouble c_even, ddouble c_odd) BesselJYCoef(ddouble x, int terms = 35) {
                    ddouble v = 1 / (x * x), w = -v;

                    ddouble c_even = ACoef(0), c_odd = ACoef(1);

                    for (int k = 1; k <= terms; k++) {
                        ddouble dc_even = w * ACoef(2 * k);
                        ddouble dc_odd = w * ACoef(2 * k + 1);

                        c_even += dc_even;
                        c_odd += dc_odd;

                        if (((long)ddouble.ILogB(c_even) - ddouble.ILogB(dc_even) >= 106L || ddouble.IsZero(dc_even)) &&
                            ((long)ddouble.ILogB(c_odd) - ddouble.ILogB(dc_odd) >= 106L || ddouble.IsZero(dc_odd))) {

                            break;
                        }

                        w *= -v;
                    }

                    return (c_even, c_odd / x);
                }

                public ddouble BesselICoef(ddouble x, int terms = 72) {
                    ddouble v = 1d / x, w = -v;

                    ddouble c = ACoef(0);

                    for (int k = 1; k <= terms; k++) {
                        ddouble dc = w * ACoef(k);

                        c += dc;

                        if ((long)ddouble.ILogB(c) - ddouble.ILogB(dc) >= 106L || ddouble.IsZero(dc)) {
                            break;
                        }

                        w *= -v;
                    }

                    return c;
                }

                public ddouble BesselKCoef(ddouble x, int terms = 52) {
                    ddouble v = 1d / x, w = v;

                    ddouble c = ACoef(0);

                    for (int k = 1; k <= terms; k++) {
                        ddouble dc = w * ACoef(k);

                        c += dc;

                        if ((long)ddouble.ILogB(c) - ddouble.ILogB(dc) >= 106L || ddouble.IsZero(dc)) {
                            break;
                        }

                        w *= v;
                    }

                    return c;
                }
            }
        }

        public class MillerBackward {
            public static readonly double MillerBwdBesselYEps = double.ScaleB(1, -30);

            private static readonly Dictionary<ddouble, BesselJPhiTable> phi_coef_table = [];
            private static readonly Dictionary<ddouble, BesselYEtaTable> eta_coef_table = [];
            private static readonly Dictionary<ddouble, BesselYXiTable> xi_coef_table = [];

            public static ddouble BesselJ(int n, ddouble x) {
                Debug.Assert(ddouble.IsPositive(x));

                int m = BesselJYIterM((double)x);

                ddouble y = BesselJKernel(n, x, m);

                return y;
            }

            public static ddouble BesselJ(ddouble nu, ddouble x) {
                Debug.Assert(ddouble.IsPositive(x));

                int m = BesselJYIterM((double)x);

                if (NearlyInteger(nu, out int n)) {
                    ddouble y = BesselJKernel(n, x, m);

                    return y;
                }
                else {
                    ddouble y = BesselJKernel(nu, x, m);

                    return y;
                }
            }

            public static ddouble BesselY(int n, ddouble x) {
                Debug.Assert(ddouble.IsPositive(x));

                int m = BesselJYIterM((double)x);

                ddouble y = BesselYKernel(n, x, m);

                return y;
            }

            public static ddouble BesselY(ddouble nu, ddouble x) {
                Debug.Assert(ddouble.IsPositive(x));

                int m = BesselJYIterM((double)x);

                if (NearlyInteger(nu, out int n)) {
                    ddouble y = BesselYKernel(n, x, m);

                    return y;
                }
                else {
                    ddouble y = BesselYKernel(nu, x, m);

                    return y;
                }
            }

            private static int BesselJYIterM(double r) {
                int m = (int)double.Ceiling(3.8029e1 + r * 1.6342e0);

                return (m + 1) / 2 * 2;
            }

            private static ddouble BesselJKernel(int n, ddouble x, int m) {
                Debug.Assert(m >= 2 && (m & 1) == 0 && n < m);

                if (n < 0) {
                    return (n & 1) == 0 ? BesselJKernel(-n, x, m) : -BesselJKernel(-n, x, m);
                }
                else if (n == 0) {
                    return BesselJ0Kernel(x, m);
                }
                else if (n == 1) {
                    return BesselJ1Kernel(x, m);
                }
                else {
                    return BesselJNKernel(n, x, m);
                }
            }

            private static ddouble BesselJKernel(ddouble nu, ddouble x, int m) {
                int n = (int)ddouble.Floor(nu);
                ddouble alpha = nu - n;

                Debug.Assert(m >= 2 && (m & 1) == 0 && n < m);

                if (!phi_coef_table.TryGetValue(alpha, out BesselJPhiTable phi)) {
                    phi = new BesselJPhiTable(alpha);
                    phi_coef_table.Add(alpha, phi);
                }

                ddouble f0 = 1e-256d, f1 = 0d, lambda = 0d;
                ddouble v = 1d / x;

                if (n >= 0) {
                    ddouble fn = 0d;

                    for (int k = m; k >= 1; k--) {
                        if ((k & 1) == 0) {
                            lambda += f0 * phi[k / 2];
                        }

                        (f0, f1) = (ddouble.Ldexp(k + alpha, 1) * v * f0 - f1, f0);

                        if (k - 1 == n) {
                            fn = f0;
                        }
                    }

                    lambda += f0 * phi[0];
                    lambda *= ddouble.Pow(ddouble.Ldexp(v, 1), alpha);

                    ddouble yn = fn / lambda;

                    return yn;
                }
                else {
                    for (int k = m; k >= 1; k--) {
                        if ((k & 1) == 0) {
                            lambda += f0 * phi[k / 2];
                        }

                        (f0, f1) = (ddouble.Ldexp(k + alpha, 1) * v * f0 - f1, f0);
                    }

                    lambda += f0 * phi[0];
                    lambda *= ddouble.Pow(ddouble.Ldexp(v, 1), alpha);

                    for (int k = 0; k > n; k--) {
                        (f0, f1) = (ddouble.Ldexp(k + alpha, 1) * v * f0 - f1, f0);
                    }

                    ddouble yn = f0 / lambda;

                    return yn;
                }
            }

            private static ddouble BesselJ0Kernel(ddouble x, int m) {
                Debug.Assert(m >= 2 && (m & 1) == 0);

                ddouble f0 = 1e-256d, f1 = 0d, lambda = 0d;
                ddouble v = 1d / x;

                for (int k = m; k >= 1; k--) {
                    if ((k & 1) == 0) {
                        lambda += f0;
                    }

                    (f0, f1) = (2 * k * v * f0 - f1, f0);
                }

                lambda = ddouble.Ldexp(lambda, 1) + f0;

                ddouble y0 = f0 / lambda;

                return y0;
            }

            private static ddouble BesselJ1Kernel(ddouble x, int m) {
                Debug.Assert(m >= 2 && (m & 1) == 0);

                ddouble f0 = 1e-256d, f1 = 0d, lambda = 0d;
                ddouble v = 1d / x;

                for (int k = m; k >= 1; k--) {
                    if ((k & 1) == 0) {
                        lambda += f0;
                    }

                    (f0, f1) = (2 * k * v * f0 - f1, f0);
                }

                lambda = ddouble.Ldexp(lambda, 1) + f0;

                ddouble y1 = f1 / lambda;

                return y1;
            }

            private static ddouble BesselJNKernel(int n, ddouble x, int m) {
                Debug.Assert(m >= 2 && (m & 1) == 0);

                ddouble f0 = 1e-256d, f1 = 0d, fn = 0d, lambda = 0d;
                ddouble v = 1d / x;

                for (int k = m; k >= 1; k--) {
                    if ((k & 1) == 0) {
                        lambda += f0;
                    }

                    (f0, f1) = (2 * k * v * f0 - f1, f0);

                    if (k - 1 == n) {
                        fn = f0;
                    }
                }

                lambda = ddouble.Ldexp(lambda, 1) + f0;

                ddouble yn = fn / lambda;

                return yn;
            }

            private static ddouble BesselYKernel(int n, ddouble x, int m) {
                Debug.Assert(m >= 2 && (m & 1) == 0 && n < m);

                if (n < 0) {
                    return (n & 1) == 0 ? BesselYKernel(-n, x, m) : -BesselYKernel(-n, x, m);
                }
                else if (n == 0) {
                    return BesselY0Kernel(x, m);
                }
                else if (n == 1) {
                    return BesselY1Kernel(x, m);
                }
                else {
                    return BesselYNKernel(n, x, m);
                }
            }

            private static ddouble BesselYKernel(ddouble nu, ddouble x, int m) {
                int n = (int)ddouble.Round(nu);
                ddouble alpha = nu - n;

                Debug.Assert(m >= 2 && (m & 1) == 0 && n < m);

                if (!eta_coef_table.TryGetValue(alpha, out BesselYEtaTable eta)) {
                    eta = new BesselYEtaTable(alpha);
                    eta_coef_table.Add(alpha, eta);
                }
                if (!xi_coef_table.TryGetValue(alpha, out BesselYXiTable xi)) {
                    xi = new BesselYXiTable(alpha, eta);
                    xi_coef_table.Add(alpha, xi);
                }
                if (!phi_coef_table.TryGetValue(alpha, out BesselJPhiTable phi)) {
                    phi = new BesselJPhiTable(alpha);
                    phi_coef_table.Add(alpha, phi);
                }

                ddouble f0 = 1e-256, f1 = 0d, lambda = 0d;
                ddouble se = 0d, sxo = 0d, sxe = 0d;
                ddouble v = 1d / x;

                for (int k = m; k >= 1; k--) {
                    if ((k & 1) == 0) {
                        lambda += f0 * phi[k / 2];

                        se += f0 * eta[k / 2];
                        sxe += f0 * xi[k];
                    }
                    else if (k >= 3) {
                        sxo += f0 * xi[k];
                    }

                    (f0, f1) = (ddouble.Ldexp(k + alpha, 1) * v * f0 - f1, f0);
                }

                ddouble s = ddouble.Pow(ddouble.Ldexp(v, 1), alpha), sqs = s * s;

                lambda += f0 * phi[0];
                lambda *= s;

                ddouble rcot = 1d / ddouble.TanPI(alpha), rgamma = ddouble.Gamma(1d + alpha), rsqgamma = rgamma * rgamma;
                ddouble r = ddouble.Ldexp(ddouble.RcpPI * sqs, 1);
                ddouble p = sqs * rsqgamma * ddouble.RcpPI;

                ddouble eta0 = ddouble.Abs(alpha) > MillerBwdBesselYEps
                    ? rcot - p / alpha
                    : BesselYEta0Eps(alpha, x);

                ddouble xi0 = -ddouble.Ldexp(v, 1) * p;
                ddouble xi1 = ddouble.Abs(alpha) > MillerBwdBesselYEps
                    ? rcot + p * (alpha * (alpha + 1d) + 1d) / (alpha * (alpha - 1d))
                    : BesselYXi1Eps(alpha, x);

                ddouble y0 = r * se + eta0 * f0;
                ddouble y1 = r * (3d * alpha * v * sxe + sxo) + xi0 * f0 + xi1 * f1;

                if (n == 0) {
                    ddouble yn = y0 / lambda;

                    return yn;
                }
                if (n == 1) {
                    ddouble yn = y1 / lambda;

                    return yn;
                }
                if (n >= 0) {
                    for (int k = 1; k < n; k++) {
                        (y1, y0) = (ddouble.Ldexp(k + alpha, 1) * v * y1 - y0, y1);
                    }

                    ddouble yn = y1 / lambda;

                    return yn;
                }
                else {
                    for (int k = 0; k > n; k--) {
                        (y0, y1) = (ddouble.Ldexp(k + alpha, 1) * v * y0 - y1, y0);
                    }

                    ddouble yn = y0 / lambda;

                    return yn;
                }
            }

            private static ddouble BesselYEta0Eps(ddouble alpha, ddouble x) {
                ddouble lnx = ddouble.Log(x), lnhalfx = ddouble.Log(ddouble.Ldexp(x, -1));
                ddouble pi = ddouble.PI, sqpi = pi * pi;
                ddouble ln2 = ddouble.Ln2, sqln2 = ln2 * ln2, cbln2 = sqln2 * ln2, qdln2 = sqln2 * sqln2;
                ddouble g = ddouble.EulerGamma;

                ddouble r0 = lnhalfx + g;
                ddouble r1 =
                    (-sqln2 + lnx * (ln2 * 2d - lnx)) * 4d
                    - sqpi
                    - g * (lnhalfx * 2d + g) * 4d;
                ddouble r2 =
                    (-cbln2 + lnx * (sqln2 * 3d + lnx * (ln2 * -3d + lnx))) * 4d
                    + ddouble.Zeta3 * 2d
                    + sqpi * (lnhalfx + g)
                    + g * ((sqln2 + lnx * (ln2 * -2d + lnx)) * 3d + g * (lnhalfx * 3d + g)) * 4d;
                ddouble r3 =
                    (-qdln2 + lnx * (cbln2 * 4d + lnx * (sqln2 * -6d + lnx * (ln2 * 4d - lnx)))) * 16d
                    - ddouble.Zeta3 * (lnhalfx + g) * 32d
                    - sqpi * ((sqln2 + lnx * (-ln2 * 2d + lnx) + g * (lnhalfx * 2d + g)) * 8d + sqpi)
                    + g * ((cbln2 + lnx * (sqln2 * -3d + lnx * (ln2 * 3d - lnx))) * 4d
                    + g * ((sqln2 + lnx * (ln2 * -2d + lnx)) * -6d
                    + g * (lnhalfx * -4d
                    - g))) * 16d;

                ddouble eta0 = (r0 * 48d + alpha * (r1 * 12d + alpha * (r2 * 8d + alpha * r3))) / (24d * ddouble.PI);

                return eta0;
            }

            private static ddouble BesselYXi1Eps(ddouble alpha, ddouble x) {
                ddouble lnx = ddouble.Log(x), lnhalfx = ddouble.Log(ddouble.Ldexp(x, -1)), lnxm1 = lnx - 1, lnhalfxm1 = lnhalfx - 1;
                ddouble pi = ddouble.PI, sqpi = pi * pi;
                ddouble ln2 = ddouble.Ln2, sqln2 = ln2 * ln2, cbln2 = sqln2 * ln2, qdln2 = sqln2 * sqln2;
                ddouble g = ddouble.EulerGamma;

                ddouble r0 = lnhalfxm1 + g;
                ddouble r1 =
                    (-sqln2 + ln2 * lnxm1 * 2d + lnx * (2 - lnx)) * 4d
                    - sqpi
                    - g * (lnhalfxm1 * 2d + g) * 4d
                    - 6d;
                ddouble r2 =
                    -cbln2 * 4d + sqln2 * lnxm1 * 12d + lnx * (18d + lnx * (-12d + lnx * 4d))
                    + ln2 * (lnx * (2d - lnx) * 12d - 18d)
                    + ddouble.Zeta3 * 2d
                    + sqpi * (lnhalfxm1 + g)
                    + g * ((sqln2 - ln2 * lnxm1 * 2d + lnx * (-2d + lnx)) * 12d + 18d
                    + g * (lnhalfxm1 * 12d
                    + g * 4d))
                    - 9d;
                ddouble r3 =
                    -qdln2 * 16d
                    + cbln2 * lnxm1 * 64d
                    + sqln2 * (lnx * (2d - lnx) * 96d - 144d)
                    + ln2 * (lnx * (9d + lnx * (-6d + lnx * 2d)) * 32d - 144d)
                    + lnx * (9d + lnx * (-9d + lnx * (4d - lnx))) * 16d
                    + ddouble.Zeta3 * (lnhalfxm1 + g) * -32d
                    + sqpi * ((-sqln2 + ln2 * lnxm1 * 2d + lnx * (2d - lnx) - g * (lnhalfxm1 * 2d + g)) * 8d - 12d - sqpi)
                    + g * ((cbln2 - sqln2 * lnxm1 * 3d) * 64d + ln2 * (lnx * (-2d + lnx) * 192d + 288d) + lnx * (-9d + lnx * (6d - lnx * 2d)) * 32d + 144d
                    + g * ((-sqln2 + ln2 * lnxm1 * 2d + lnx * (2d - lnx)) * 96d - 144d
                    + g * (lnhalfxm1 * -64d
                    - g * 16d)))
                    - 72d;

                ddouble xi1 = (r0 * 48d + alpha * (r1 * 12d + alpha * (r2 * 8d + alpha * r3))) / (24d * ddouble.PI);

                return xi1;
            }

            private static ddouble BesselY0Kernel(ddouble x, int m) {
                Debug.Assert(m >= 2 && (m & 1) == 0);

                if (!eta_coef_table.TryGetValue(0, out BesselYEtaTable eta)) {
                    eta = new BesselYEtaTable(0);
                    eta_coef_table.Add(0, eta);
                }

                ddouble f0 = 1e-256, f1 = 0d, lambda = 0d;
                ddouble se = 0d;
                ddouble v = 1d / x;

                for (int k = m; k >= 1; k--) {
                    if ((k & 1) == 0) {
                        lambda += f0;

                        se += f0 * eta[k / 2];
                    }

                    (f0, f1) = (2 * k * v * f0 - f1, f0);
                }

                lambda = ddouble.Ldexp(lambda, 1) + f0;

                ddouble y0 = ddouble.Ldexp((se + f0 * (ddouble.Log(ddouble.Ldexp(x, -1)) + ddouble.EulerGamma)) / (ddouble.PI * lambda), 1);

                return y0;
            }

            private static ddouble BesselY1Kernel(ddouble x, int m) {
                Debug.Assert(m >= 2 && (m & 1) == 0);

                if (!xi_coef_table.TryGetValue(0, out BesselYXiTable xi)) {
                    if (!eta_coef_table.TryGetValue(0, out BesselYEtaTable eta)) {
                        eta = new BesselYEtaTable(0);
                        eta_coef_table.Add(0, eta);
                    }

                    xi = new BesselYXiTable(0, eta);
                    xi_coef_table.Add(0, xi);
                }

                ddouble f0 = 1e-256, f1 = 0d, lambda = 0d;
                ddouble sx = 0d;
                ddouble v = 1d / x;

                for (int k = m; k >= 1; k--) {
                    if ((k & 1) == 0) {
                        lambda += f0;
                    }
                    else if (k >= 3) {
                        sx += f0 * xi[k];
                    }

                    (f0, f1) = (2 * k * v * f0 - f1, f0);
                }

                lambda = ddouble.Ldexp(lambda, 1) + f0;

                ddouble y1 = ddouble.Ldexp((sx - v * f0 + (ddouble.Log(ddouble.Ldexp(x, -1)) + ddouble.EulerGamma - 1d) * f1) / (lambda * ddouble.PI), 1);

                return y1;
            }

            private static ddouble BesselYNKernel(int n, ddouble x, int m) {
                Debug.Assert(m >= 2 && (m & 1) == 0);

                if (!eta_coef_table.TryGetValue(0, out BesselYEtaTable eta)) {
                    eta = new BesselYEtaTable(0);
                    eta_coef_table.Add(0, eta);
                }
                if (!xi_coef_table.TryGetValue(0, out BesselYXiTable xi)) {
                    xi = new BesselYXiTable(0, eta);
                    xi_coef_table.Add(0, xi);
                }

                ddouble f0 = 1e-256, f1 = 0d, lambda = 0d;
                ddouble se = 0d, sx = 0d;
                ddouble v = 1d / x;

                for (int k = m; k >= 1; k--) {
                    if ((k & 1) == 0) {
                        lambda += f0;

                        se += f0 * eta[k / 2];
                    }
                    else if (k >= 3) {
                        sx += f0 * xi[k];
                    }

                    (f0, f1) = (2 * k * v * f0 - f1, f0);
                }

                lambda = ddouble.Ldexp(lambda, 1) + f0;

                ddouble c = ddouble.Log(ddouble.Ldexp(x, -1)) + ddouble.EulerGamma;

                ddouble y0 = se + f0 * c;
                ddouble y1 = sx - v * f0 + (c - 1d) * f1;

                for (int k = 1; k < n; k++) {
                    (y1, y0) = (2 * k * v * y1 - y0, y1);
                }

                ddouble yn = ddouble.Ldexp(y1 / (lambda * ddouble.PI), 1);

                return yn;
            }

            private class BesselJPhiTable {
                private readonly ddouble alpha;
                private readonly List<ddouble> table = [];

                private ddouble g;

                public BesselJPhiTable(ddouble alpha) {
                    Debug.Assert(alpha > -1d && alpha < 1d, nameof(alpha));

                    this.alpha = alpha;

                    ddouble phi0 = ddouble.Gamma(1 + alpha);
                    ddouble phi1 = phi0 * (alpha + 2d);

                    this.g = phi0;

                    this.table.Add(phi0);
                    this.table.Add(phi1);
                }

                public ddouble this[int k] => Value(k);

                private ddouble Value(int k) {
                    Debug.Assert(k >= 0);

                    if (k < table.Count) {
                        return table[k];
                    }

                    for (int i = table.Count; i <= k; i++) {
                        g = g * (alpha + i - 1d) / i;

                        ddouble phi = g * (alpha + 2 * i);

                        table.Add(phi);
                    }

                    return table[k];
                }
            };

            private class BesselIPsiTable {
                private readonly ddouble alpha;
                private readonly List<ddouble> table = [];

                private ddouble g;

                public BesselIPsiTable(ddouble alpha) {
                    Debug.Assert(alpha > -1d && alpha < 1d, nameof(alpha));

                    this.alpha = alpha;

                    ddouble psi0 = ddouble.Gamma(1d + alpha);
                    ddouble psi1 = ddouble.Ldexp(psi0, 1) * (1d + alpha);

                    this.g = ddouble.Ldexp(psi0, 1);

                    this.table.Add(psi0);
                    this.table.Add(psi1);
                }

                public ddouble this[int k] => Value(k);

                private ddouble Value(int k) {
                    Debug.Assert(k >= 0);

                    if (k < table.Count) {
                        return table[k];
                    }

                    for (int i = table.Count; i <= k; i++) {
                        g = g * (ddouble.Ldexp(alpha, 1) + i - 1d) / i;

                        ddouble phi = g * (alpha + i);

                        table.Add(phi);
                    }

                    return table[k];
                }
            };

            private class BesselYEtaTable {
                private readonly ddouble alpha;
                private readonly List<ddouble> table = [];

                private ddouble g;

                public BesselYEtaTable(ddouble alpha) {
                    Debug.Assert(alpha > -1d && alpha < 1d, nameof(alpha));

                    this.alpha = alpha;
                    this.table.Add(ddouble.NaN);

                    if (alpha != 0d) {
                        ddouble c = ddouble.Gamma(1d + alpha);
                        c *= c;
                        this.g = 1d / (1d - alpha) * c;

                        ddouble eta1 = (alpha + 2d) * g;

                        this.table.Add(eta1);
                    }
                }

                public ddouble this[int k] => Value(k);

                private ddouble Value(int k) {
                    Debug.Assert(k >= 0);

                    if (k < table.Count) {
                        return table[k];
                    }

                    for (int i = table.Count; i <= k; i++) {
                        if (alpha != 0d) {
                            g = -g * (alpha + i - 1) * (ddouble.Ldexp(alpha, 1) + i - 1d) / (i * (i - alpha));

                            ddouble eta = g * (alpha + 2 * i);

                            table.Add(eta);
                        }
                        else {
                            ddouble eta = (ddouble)2d / i;

                            table.Add((i & 1) == 1 ? eta : -eta);
                        }
                    }

                    return table[k];
                }
            };

            private class BesselYXiTable {
                private readonly ddouble alpha;
                private readonly List<ddouble> table = [];
                private readonly BesselYEtaTable eta;

                public BesselYXiTable(ddouble alpha, BesselYEtaTable eta) {
                    Debug.Assert(alpha >= -1d && alpha < 1d, nameof(alpha));

                    this.alpha = alpha;
                    this.table.Add(ddouble.NaN);
                    this.table.Add(ddouble.NaN);

                    this.eta = eta;
                }

                public ddouble this[int k] => Value(k);

                private ddouble Value(int k) {
                    Debug.Assert(k >= 0);

                    if (k < table.Count) {
                        return table[k];
                    }

                    for (int i = table.Count; i <= k; i++) {
                        if (alpha != 0d) {
                            if ((i & 1) == 0) {
                                table.Add(eta[i / 2]);
                            }
                            else {
                                table.Add((eta[i / 2] - eta[i / 2 + 1]) / 2);
                            }
                        }
                        else {
                            if ((i & 1) == 1) {
                                ddouble xi = (ddouble)(2 * (i / 2) + 1) / (i / 2 * (i / 2 + 1));
                                table.Add((i & 2) > 0 ? xi : -xi);
                            }
                            else {
                                table.Add(ddouble.NaN);
                            }
                        }
                    }

                    return table[k];
                }
            }
        }

        public static class YoshidaPade {
            private static readonly ReadOnlyCollection<ReadOnlyCollection<ddouble>> ess_coef_table;
            private static readonly Dictionary<ddouble, ReadOnlyCollection<(ddouble c, ddouble s)>> cds_coef_table = [];

            static YoshidaPade() {
                cds_coef_table.Add(0, Array.AsReadOnly(YoshidaPadeCoefM31.Nu0.Reverse().ToArray()));
                cds_coef_table.Add(1, Array.AsReadOnly(YoshidaPadeCoefM31.Nu1.Reverse().ToArray()));

                List<ReadOnlyCollection<ddouble>> es = [];

                for (int i = 0; i < YoshidaPadeCoefM31.Ess.Length; i++) {
                    es.Add(Array.AsReadOnly(YoshidaPadeCoefM31.Ess[i].ToArray()));
                }

                ess_coef_table = Array.AsReadOnly(es.ToArray());
            }

            public static ddouble BesselK(ddouble nu, ddouble x, bool scale = false) {
                if (nu < 2d) {
                    if (!cds_coef_table.TryGetValue(nu, out ReadOnlyCollection<(ddouble c, ddouble s)> cds_table)) {
                        cds_table = Table(nu);
                        cds_coef_table.Add(nu, cds_table);
                    }

                    ReadOnlyCollection<(ddouble, ddouble)> cds = cds_table;

                    ddouble y = Value(x, cds, scale);

                    return y;
                }
                else {
                    int n = (int)ddouble.Floor(nu);
                    ddouble alpha = nu - n;

                    ddouble y0 = BesselK(alpha, x, scale);
                    ddouble y1 = BesselK(alpha + 1d, x, scale);

                    ddouble v = 1d / x;

                    for (int k = 1; k < n; k++) {
                        (y1, y0) = (ddouble.Ldexp(k + alpha, 1) * v * y1 + y0, y1);
                    }

                    return y1;
                }
            }

            private static ddouble Value(ddouble x, ReadOnlyCollection<(ddouble c, ddouble d)> cds, bool scale = false) {
                ddouble t = 1d / x;
                (ddouble sc, ddouble sd) = cds[0];

                for (int i = 1; i < cds.Count; i++) {
                    (ddouble c, ddouble d) = cds[i];

                    sc = sc * t + c;
                    sd = sd * t + d;
                }

                ddouble y = ddouble.Sqrt(ddouble.Ldexp(t * ddouble.PI, -1)) * sc / sd;

                if (!scale) {
                    y *= ddouble.Exp(-x);
                }

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

        public static class CubicInterpolate {
            public static ddouble BesselYPowerSeries(ddouble nu, ddouble x) {
                int n = (int)ddouble.Round(nu);
                ddouble alpha = nu - n;

                ddouble y0 = PowerSeries.BesselY(n, x);
                ddouble y1 = PowerSeries.BesselY(n + ddouble.Sign(alpha) * InterpolationThreshold, x);
                ddouble y2 = PowerSeries.BesselY(n + ddouble.Sign(alpha) * InterpolationThreshold * 1.5d, x);
                ddouble y3 = PowerSeries.BesselY(n + ddouble.Sign(alpha) * InterpolationThreshold * 2d, x);

                ddouble t = ddouble.Abs(alpha) / InterpolationThreshold;
                ddouble y = Interpolate(t, y0, y1, y2, y3);

                return y;
            }

            public static ddouble BesselYMillerBackward(ddouble nu, ddouble x) {
                int n = (int)ddouble.Round(nu);
                ddouble alpha = nu - n;

                ddouble y0 = MillerBackward.BesselY(n, x);
                ddouble y1 = MillerBackward.BesselY(n + ddouble.Sign(alpha) * InterpolationThreshold, x);
                ddouble y2 = MillerBackward.BesselY(n + ddouble.Sign(alpha) * InterpolationThreshold * 1.5d, x);
                ddouble y3 = MillerBackward.BesselY(n + ddouble.Sign(alpha) * InterpolationThreshold * 2d, x);

                ddouble t = ddouble.Abs(alpha) / InterpolationThreshold;
                ddouble y = Interpolate(t, y0, y1, y2, y3);

                return y;
            }

            public static ddouble BesselKPowerSeries(ddouble nu, ddouble x, bool scale = false) {
                int n = (int)ddouble.Round(nu);
                ddouble alpha = nu - n;

                ddouble y0 = PowerSeries.BesselK(n, x, scale);
                ddouble y1 = PowerSeries.BesselK(n + ddouble.Sign(alpha) * InterpolationThreshold, x, scale);
                ddouble y2 = PowerSeries.BesselK(n + ddouble.Sign(alpha) * InterpolationThreshold * 1.5d, x, scale);
                ddouble y3 = PowerSeries.BesselK(n + ddouble.Sign(alpha) * InterpolationThreshold * 2d, x, scale);

                ddouble t = ddouble.Abs(alpha) / InterpolationThreshold;
                ddouble y = Interpolate(t, y0, y1, y2, y3);

                return y;
            }

            private static ddouble Interpolate(ddouble t, ddouble y0, ddouble y1, ddouble y2, ddouble y3) {
                return y0 + (
                    -(13d + t * (-9d + t * 2d)) / 6d * y0
                    + (6d + t * (-7d + t * 2d)) * y1
                    - (16d + t * (-24d + t * 8d)) / 3d * y2
                    + (3d + t * (-5d + t * 2d)) / 2d * y3) * t;
            }
        }
    }
}