using DoubleDouble;
using DoubleDoubleComplex;
using System.Collections.ObjectModel;
using System.Diagnostics;
using static DDoubleOptimizedBessel.ComplexBesselUtil;

namespace DDoubleOptimizedBessel {

    public static class ComplexBessel {

        public static Complex BesselJ(ddouble nu, Complex z) {
            CheckNu(nu);

            if (ddouble.IsNegative(z.I)) {
                return BesselJ(nu, z.Conj).Conj;
            }

            if (ddouble.IsNegative(z.R)) {
                return (SinCosPICache.CosPI(nu), SinCosPICache.SinPI(nu)) * BesselJ(nu, -z);
            }

            if (z.Magnitude >= HankelThreshold) {
                return Limit.BesselJ(nu, z);
            }
            else if (z.R <= PowerSeriesThreshold(nu, z.I)) {
                return PowerSeries.BesselJ(nu, z);
            }
            else if (z.I <= MillerBackwardThreshold) {
                return MillerBackward.BesselJ(nu, z);
            }
            else {
                return (SinCosPICache.CosPI(nu / 2), SinCosPICache.SinPI(nu / 2)) * MillerBackward.BesselI(nu, (z.I, z.R)).Conj;
            }
        }

        public static Complex BesselY(ddouble nu, Complex z) {
            CheckNu(nu);

            if (ddouble.IsNegative(z.I)) {
                return BesselY(nu, z.Conj).Conj;
            }

            if (ddouble.IsNegative(z.R)) {
                return (SinCosPICache.CosPI(nu), -SinCosPICache.SinPI(nu)) * BesselY(nu, -z)
                        + (0, 2 * SinCosPICache.CosPI(nu)) * BesselJ(nu, -z);
            }

            if (z.Magnitude >= HankelThreshold) {
                return Limit.BesselY(nu, z);
            }
            else if (z.R <= PowerSeriesThreshold(nu, z.I) - BesselJYPowerseriesBias) {
                if (NearlyInteger(nu, out _) || ddouble.Abs(ddouble.Round(nu) - nu) >= InterpolationThreshold) {
                    return PowerSeries.BesselY(nu, z);
                }
                else {
                    return CubicInterpolate.BesselYPowerSeries(nu, z);
                }
            }
            else if (z.I <= MillerBackwardThreshold) {
                if (NearlyInteger(nu, out _) || ddouble.Abs(ddouble.Ceiling(nu) - nu) >= InterpolationThreshold) {
                    return MillerBackward.BesselY(nu, z);
                }
                else {
                    return CubicInterpolate.BesselYMillerBackward(nu, z);
                }
            }
            else {
                Complex c = (SinCosPICache.CosPI(nu / 2), SinCosPICache.SinPI(nu / 2));

                Complex bi = (z.R <= PowerSeriesThreshold(nu, z.I))
                    ? PowerSeries.BesselI(nu, (z.I, z.R))
                    : MillerBackward.BesselI(nu, (z.I, z.R));

                Complex bk = YoshidaPade.BesselK(nu, (z.I, z.R));

                Complex y = (0, 1) * c * bi.Conj - 2 * ddouble.RcpPI * (c * bk).Conj;

                return y;
            }
        }

        public static Complex BesselI(ddouble nu, Complex z) {
            CheckNu(nu);

            if (ddouble.IsNegative(z.I)) {
                return BesselI(nu, z.Conj).Conj;
            }

            if (ddouble.IsNegative(z.R)) {
                return (SinCosPICache.CosPI(nu), SinCosPICache.SinPI(nu)) * BesselI(nu, -z);
            }

            if (z.Magnitude >= HankelThreshold) {
                return Limit.BesselI(nu, z);
            }
            else if (z.I <= PowerSeriesThreshold(nu, z.R)) {
                return PowerSeries.BesselI(nu, z);
            }
            else {
                return MillerBackward.BesselI(nu, z);
            }
        }

        public static Complex BesselK(ddouble nu, Complex z) {
            CheckNu(nu);

            nu = ddouble.Abs(nu);

            if (ddouble.IsNegative(z.I)) {
                return BesselK(nu, z.Conj).Conj;
            }

            if (ddouble.IsNegative(z.R)) {
                return (SinCosPICache.CosPI(nu), -SinCosPICache.SinPI(nu)) * BesselK(nu, -z)
                        - (0, ddouble.PI) * BesselI(nu, -z);
            }

            if (z.Magnitude >= HankelThreshold) {
                return Limit.BesselK(nu, z);
            }
            else if (z.Magnitude <= BesselKNearZeroThreshold) {
                if (NearlyInteger(nu, out _) || ddouble.Abs(ddouble.Round(nu) - nu) >= InterpolationThreshold) {
                    return PowerSeries.BesselK(nu, z);
                }
                else {
                    return CubicInterpolate.BesselKPowerSeries(nu, z);
                }
            }
            else if (z.R >= BesselKPadeThreshold) {
                return YoshidaPade.BesselK(nu, z);
            }
            else {
                Complex c = (SinCosPICache.CosPI(nu / 2), -SinCosPICache.SinPI(nu / 2));

                Complex bi = (z.I <= PowerSeriesThreshold(nu, z.R))
                    ? PowerSeries.BesselI(nu, z)
                    : MillerBackward.BesselI(nu, z);

                Complex by = (z.I <= PowerSeriesThreshold(nu, z.R) - BesselJYPowerseriesBias)
                    ? ((NearlyInteger(nu, out _) || ddouble.Abs(ddouble.Round(nu) - nu) >= InterpolationThreshold)
                        ? PowerSeries.BesselY(nu, (z.I, z.R)) : CubicInterpolate.BesselYPowerSeries(nu, (z.I, z.R))
                    )
                    : ((NearlyInteger(nu, out _) || ddouble.Abs(ddouble.Ceiling(nu) - nu) >= InterpolationThreshold)
                        ? MillerBackward.BesselY(nu, (z.I, z.R)) : CubicInterpolate.BesselYMillerBackward(nu, (z.I, z.R))
                    );

                Complex y = c * ((0, -1) * c * bi - by.Conj) * ddouble.PI / 2;

                return y;
            }
        }
    }

    public static class ComplexBesselUtil {
        public const int MaxN = 16;
        public static readonly double Eps = double.ScaleB(1, -1000);
        public static readonly double ExtremelyNearZero = double.ScaleB(1, -28);
        public static readonly double InterpolationThreshold = double.ScaleB(1, -25);
        public const double HankelThreshold = 38.75, MillerBackwardThreshold = 6;
        public const double BesselKPadeThreshold = 1, BesselKNearZeroThreshold = 4, BesselJYPowerseriesBias = 2;

        public static ddouble PowerSeriesThreshold(ddouble nu, ddouble x) {
            ddouble nu_abs = ddouble.Abs(nu);

            return 7.5 + nu_abs * (3.57e-1 + nu_abs * 5.23e-3) + x * (4.67e-1 - nu_abs * 1.51e-2);
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

            public static Complex BesselJ(ddouble nu, Complex z) {
                Debug.Assert(ddouble.IsPositive(z.R));
                Debug.Assert(ddouble.IsPositive(z.I));

                if (ddouble.IsNegative(nu) && NearlyInteger(nu, out int n)) {
                    Complex y = BesselJ(-nu, z);

                    return (n & 1) == 0 ? y : -y;
                }
                else {
                    Complex y = BesselJIKernel(nu, z, sign_switch: true, terms: 42);

                    return y;
                }
            }

            public static Complex BesselY(ddouble nu, Complex z) {
                Debug.Assert(ddouble.IsPositive(z.R));
                Debug.Assert(ddouble.IsPositive(z.I));

                if (NearlyInteger(nu, out int n)) {
                    Complex y = BesselYKernel(n, z, terms: 44);

                    return y;
                }
                else {
                    Complex y = BesselYKernel(nu, z, terms: 44);

                    return y;
                }
            }

            public static Complex BesselI(ddouble nu, Complex z, bool scale = false) {
                Debug.Assert(ddouble.IsPositive(z.R));
                Debug.Assert(ddouble.IsPositive(z.I));

                if (ddouble.IsNegative(nu) && NearlyInteger(nu, out _)) {
                    Complex y = BesselI(-nu, z);

                    if (scale) {
                        y *= Complex.Exp(-z);
                    }

                    return y;
                }
                else {
                    Complex y = BesselJIKernel(nu, z, sign_switch: false, terms: 42);

                    if (scale) {
                        y *= Complex.Exp(-z);
                    }

                    return y;
                }
            }

            public static Complex BesselK(ddouble nu, Complex z, bool scale = false) {
                Debug.Assert(ddouble.IsPositive(z.R));
                Debug.Assert(ddouble.IsPositive(z.I));

                if (NearlyInteger(nu, out int n)) {
                    Complex y = BesselKKernel(n, z, terms: 27);

                    if (scale) {
                        y *= Complex.Exp(z);
                    }

                    return y;
                }
                else {
                    Complex y = BesselKKernel(nu, z, terms: 30);

                    if (scale) {
                        y *= Complex.Exp(z);
                    }

                    return y;
                }
            }

            private static Complex BesselJIKernel(ddouble nu, Complex z, bool sign_switch, int terms) {
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

                for (int k = 0, conv_times = 0; k <= terms && conv_times < 2; k++) {
                    Complex w = z2 * d[k];
                    Complex dc = u * r[k] * (sign_switch ? 1d - w : 1d + w);

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

                for (int k = 0, t = 1, conv_times = 0; k <= terms && conv_times < 2; k++, t += 2) {
                    Complex a = t * s * g[t], q = gpn[t];
                    Complex pa = p / a, qa = q / a;

                    Complex dc = u * r[k] * (4 * t * nu * (pa + qa) - (z2 - 4 * t * t) * (pa - qa));

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
                    this.r = ddouble.Pow(4, nu);
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

            public static Complex BesselJ(ddouble nu, Complex z) {
                Debug.Assert(ddouble.IsPositive(z.R));
                Debug.Assert(ddouble.IsPositive(z.I));

                if (!table.TryGetValue(nu, out HankelExpansion hankel)) {
                    hankel = new HankelExpansion(nu);
                    table.Add(nu, hankel);
                }

                (Complex c_even, Complex c_odd) = hankel.BesselJYCoef(z);

                Complex omega = hankel.Omega(z);

                Complex cos = Complex.Cos(omega), sin = Complex.Sin(omega);

                Complex y = Complex.Sqrt(2 / (ddouble.PI * z)) * (cos * c_even - sin * c_odd);

                return y;
            }

            public static Complex BesselY(ddouble nu, Complex z) {
                Debug.Assert(ddouble.IsPositive(z.R));
                Debug.Assert(ddouble.IsPositive(z.I));

                if (!table.TryGetValue(nu, out HankelExpansion hankel)) {
                    hankel = new HankelExpansion(nu);
                    table.Add(nu, hankel);
                }


                (Complex c_even, Complex c_odd) = hankel.BesselJYCoef(z);

                Complex omega = hankel.Omega(z);

                Complex cos = Complex.Cos(omega), sin = Complex.Sin(omega);

                Complex y = Complex.Sqrt(2 / (ddouble.PI * z)) * (sin * c_even + cos * c_odd);

                return y;
            }

            public static Complex BesselI(ddouble nu, Complex z) {
                Debug.Assert(ddouble.IsPositive(z.R));
                Debug.Assert(ddouble.IsPositive(z.I));

                if (!table.TryGetValue(nu, out HankelExpansion hankel)) {
                    hankel = new HankelExpansion(nu);
                    table.Add(nu, hankel);
                }

                Complex ci = hankel.BesselICoef(z), ck = hankel.BesselKCoef(z);

                Complex y = Complex.Sqrt(1 / (2 * ddouble.PI * z)) * (
                    Complex.Exp(z) * ci -
                    (SinCosPICache.SinPI(nu), -SinCosPICache.CosPI(nu)) * Complex.Exp(-z) * ck
                );

                return y;
            }

            public static Complex BesselK(ddouble nu, Complex z) {
                Debug.Assert(ddouble.IsPositive(nu));
                Debug.Assert(ddouble.IsPositive(z.R));
                Debug.Assert(ddouble.IsPositive(z.I));

                if (!table.TryGetValue(nu, out HankelExpansion hankel)) {
                    hankel = new HankelExpansion(nu);
                    table.Add(nu, hankel);
                }

                Complex c = hankel.BesselKCoef(z);

                Complex y = Complex.Sqrt(ddouble.PI / (2 * z)) * Complex.Exp(-z) * c;

                return y;
            }

            public class HankelExpansion {
                public ddouble Nu { get; }

                private readonly List<ddouble> a_coef;

                public HankelExpansion(ddouble nu) {
                    Nu = nu;
                    a_coef = [1];
                }

                private ddouble ACoef(int n) {
                    for (int k = a_coef.Count; k <= n; k++) {
                        ddouble a = a_coef.Last() * (4 * Nu * Nu - checked((2 * k - 1) * (2 * k - 1))) / (k * 8);
                        a_coef.Add(a);
                    }

                    return a_coef[n];
                }

                public Complex Omega(Complex z) {
                    Complex omega = z - ddouble.Ldexp(2 * Nu + 1, -2) * ddouble.PI;

                    return omega;
                }

                public (Complex c_even, Complex c_odd) BesselJYCoef(Complex z, int max_term = 35) {
                    Complex v = 1 / (z * z), w = -v;

                    Complex c_even = ACoef(0), c_odd = ACoef(1);

                    for (int k = 1; k <= max_term; k++) {
                        Complex dc_even = w * ACoef(2 * k);
                        Complex dc_odd = w * ACoef(2 * k + 1);

                        c_even += dc_even;
                        c_odd += dc_odd;

                        if (((long)Complex.ILogB(c_even) - Complex.ILogB(dc_even) >= 106L || Complex.IsZero(dc_even)) &&
                            ((long)Complex.ILogB(c_odd) - Complex.ILogB(dc_odd) >= 106L || Complex.IsZero(dc_odd))) {

                            break;
                        }

                        w *= -v;
                    }

                    return (c_even, c_odd / z);
                }

                public Complex BesselICoef(Complex z, int max_term = 81) {
                    Complex v = 1 / z, w = -v;

                    Complex c = ACoef(0);

                    for (int k = 1; k <= max_term; k++) {
                        Complex dc = w * ACoef(k);

                        c += dc;

                        if ((long)Complex.ILogB(c) - Complex.ILogB(dc) >= 106L || Complex.IsZero(dc)) {
                            break;
                        }

                        w *= -v;
                    }

                    return c;
                }

                public Complex BesselKCoef(Complex z, int max_term = 59) {
                    Complex v = 1 / z, w = v;

                    Complex c = ACoef(0);

                    for (int k = 1; k <= max_term; k++) {
                        Complex dc = w * ACoef(k);

                        c += dc;

                        if ((long)Complex.ILogB(c) - Complex.ILogB(dc) >= 106L || Complex.IsZero(dc)) {
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
            private static readonly Dictionary<ddouble, BesselIPsiTable> psi_coef_table = [];
            private static readonly Dictionary<ddouble, BesselYEtaTable> eta_coef_table = [];
            private static readonly Dictionary<ddouble, BesselYXiTable> xi_coef_table = [];

            public static Complex BesselJ(int n, Complex z) {
                Debug.Assert(ddouble.IsPositive(z.R));
                Debug.Assert(ddouble.IsPositive(z.I));

                int m = BesselJYIterM((double)z.R);

                Complex y = BesselJKernel(n, z, m);

                return y;
            }

            public static Complex BesselJ(ddouble nu, Complex z) {
                Debug.Assert(ddouble.IsPositive(z.R));
                Debug.Assert(ddouble.IsPositive(z.I));

                int m = BesselJYIterM((double)z.R);

                if (NearlyInteger(nu, out int n)) {
                    Complex y = BesselJKernel(n, z, m);

                    return y;
                }
                else {
                    Complex y = BesselJKernel(nu, z, m);

                    return y;
                }
            }

            public static Complex BesselY(int n, Complex z) {
                Debug.Assert(ddouble.IsPositive(z.R));
                Debug.Assert(ddouble.IsPositive(z.I));

                int m = BesselJYIterM((double)z.R);

                Complex y = BesselYKernel(n, z, m);

                return y;
            }

            public static Complex BesselY(ddouble nu, Complex z) {
                Debug.Assert(ddouble.IsPositive(z.R));
                Debug.Assert(ddouble.IsPositive(z.I));

                int m = BesselJYIterM((double)z.R);

                if (NearlyInteger(nu, out int n)) {
                    Complex y = BesselYKernel(n, z, m);

                    return y;
                }
                else {
                    Complex y = BesselYKernel(nu, z, m);

                    return y;
                }
            }

            private static int BesselJYIterM(double r) {
                int m = (int)double.Ceiling(3.8029e1 + r * 1.6342e0);

                return (m + 1) / 2 * 2;
            }

            public static Complex BesselI(int n, Complex z, bool scale = false) {
                Debug.Assert(ddouble.IsPositive(z.R));
                Debug.Assert(ddouble.IsPositive(z.I));

                int m = BesselIIterM((double)z.R, (double)z.I);

                Complex y = BesselIKernel(n, z, m, scale);

                return y;
            }

            public static Complex BesselI(ddouble nu, Complex z, bool scale = false) {
                Debug.Assert(ddouble.IsPositive(z.R));
                Debug.Assert(ddouble.IsPositive(z.I));

                int m = BesselIIterM((double)z.R, (double)z.I);

                if (NearlyInteger(nu, out int n)) {
                    Complex y = BesselIKernel(n, z, m, scale);

                    return y;
                }
                else {
                    Complex y = BesselIKernel(nu, z, m, scale);

                    return y;
                }
            }

            private static int BesselIIterM(double r, double i) {
                int m = (int)double.Ceiling(3.3612e1 + r * 1.3557e0 + i * 1.8485e0 - r * i * 4.3649e-2);

                return (m + 1) / 2 * 2;
            }

            private static Complex BesselJKernel(int n, Complex z, int m) {
                Debug.Assert(m >= 2 && (m & 1) == 0 && n < m);

                if (n < 0) {
                    return (n & 1) == 0 ? BesselJKernel(-n, z, m) : -BesselJKernel(-n, z, m);
                }
                if (n == 0) {
                    return BesselJ0Kernel(z, m);
                }
                if (n == 1) {
                    return BesselJ1Kernel(z, m);
                }

                Complex f0 = 1e-256d, f1 = 0d, fn = 0d, lambda = 0d;
                Complex v = 1d / z;

                for (int k = m; k >= 1; k--) {
                    if ((k & 1) == 0) {
                        lambda += f0;
                    }

                    (f0, f1) = (2 * k * v * f0 - f1, f0);

                    if (k - 1 == n) {
                        fn = f0;
                    }
                }

                lambda = Complex.Ldexp(lambda, 1) + f0;

                Complex yn = fn / lambda;

                return yn;
            }

            public static Complex BesselJKernel(ddouble nu, Complex z, int m) {
                int n = (int)ddouble.Floor(nu);
                ddouble alpha = nu - n;

                if (alpha == 0d) {
                    return BesselJKernel(n, z, m);
                }

                Debug.Assert(m >= 2 && (m & 1) == 0 && n < m);

                if (!phi_coef_table.TryGetValue(alpha, out BesselJPhiTable phi_table)) {
                    phi_table = new BesselJPhiTable(alpha);
                    phi_coef_table.Add(alpha, phi_table);
                }

                BesselJPhiTable phi = phi_table;

                Complex f0 = 1e-256d, f1 = 0d, lambda = 0d;
                Complex v = 1d / z;

                if (n >= 0) {
                    Complex fn = 0d;

                    for (int k = m; k >= 1; k--) {
                        if ((k & 1) == 0) {
                            lambda += f0 * phi[k / 2];
                        }

                        (f0, f1) = (Complex.Ldexp(k + alpha, 1) * v * f0 - f1, f0);

                        if (k - 1 == n) {
                            fn = f0;
                        }
                    }

                    lambda += f0 * phi[0];
                    lambda *= Complex.Pow(Complex.Ldexp(v, 1), alpha);

                    Complex yn = fn / lambda;

                    return yn;
                }
                else {
                    for (int k = m; k >= 1; k--) {
                        if ((k & 1) == 0) {
                            lambda += f0 * phi[k / 2];
                        }

                        (f0, f1) = (Complex.Ldexp(k + alpha, 1) * v * f0 - f1, f0);
                    }

                    lambda += f0 * phi[0];
                    lambda *= Complex.Pow(Complex.Ldexp(v, 1), alpha);

                    for (int k = 0; k > n; k--) {
                        (f0, f1) = (Complex.Ldexp(k + alpha, 1) * v * f0 - f1, f0);
                    }

                    Complex yn = f0 / lambda;

                    return yn;
                }
            }

            private static Complex BesselJ0Kernel(Complex z, int m) {
                Debug.Assert(m >= 2 && (m & 1) == 0);

                Complex f0 = 1e-256d, f1 = 0d, lambda = 0d;
                Complex v = 1d / z;

                for (int k = m; k >= 1; k--) {
                    if ((k & 1) == 0) {
                        lambda += f0;
                    }

                    (f0, f1) = (2 * k * v * f0 - f1, f0);
                }

                lambda = Complex.Ldexp(lambda, 1) + f0;

                Complex y0 = f0 / lambda;

                return y0;
            }

            private static Complex BesselJ1Kernel(Complex z, int m) {
                Debug.Assert(m >= 2 && (m & 1) == 0);

                Complex f0 = 1e-256d, f1 = 0d, lambda = 0d;
                Complex v = 1d / z;

                for (int k = m; k >= 1; k--) {
                    if ((k & 1) == 0) {
                        lambda += f0;
                    }

                    (f0, f1) = (2 * k * v * f0 - f1, f0);
                }

                lambda = Complex.Ldexp(lambda, 1) + f0;

                Complex y1 = f1 / lambda;

                return y1;
            }

            private static Complex BesselYKernel(int n, Complex z, int m) {
                Debug.Assert(m >= 2 && (m & 1) == 0 && n < m);

                if (n < 0) {
                    return (n & 1) == 0 ? BesselYKernel(-n, z, m) : -BesselYKernel(-n, z, m);
                }
                if (n == 0) {
                    return BesselY0Kernel(z, m);
                }
                if (n == 1) {
                    return BesselY1Kernel(z, m);
                }

                if (!eta_coef_table.ContainsKey(0)) {
                    eta_coef_table.Add(0, new BesselYEtaTable(0));
                }

                BesselYEtaTable eta = eta_coef_table[0];

                if (!xi_coef_table.ContainsKey(0)) {
                    xi_coef_table.Add(0, new BesselYXiTable(0, eta));
                }

                BesselYXiTable xi = xi_coef_table[0];

                Complex f0 = 1e-256, f1 = 0d, lambda = 0d;
                Complex se = 0d, sx = 0d;
                Complex v = 1d / z;

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

                lambda = Complex.Ldexp(lambda, 1) + f0;

                Complex c = Complex.Log(Complex.Ldexp(z, -1)) + ddouble.EulerGamma;

                Complex y0 = se + f0 * c;
                Complex y1 = sx - v * f0 + (c - 1d) * f1;

                for (int k = 1; k < n; k++) {
                    (y1, y0) = (2 * k * v * y1 - y0, y1);
                }

                Complex yn = Complex.Ldexp(y1 / (lambda * ddouble.PI), 1);

                return yn;
            }

            public static Complex BesselYKernel(ddouble nu, Complex z, int m) {
                int n = (int)ddouble.Floor(nu);
                ddouble alpha = nu - n;

                if (alpha == 0d) {
                    return BesselYKernel(n, z, m);
                }

                Debug.Assert(m >= 2 && (m & 1) == 0 && n < m);

                if (!eta_coef_table.TryGetValue(alpha, out BesselYEtaTable eta_table)) {
                    eta_table = new BesselYEtaTable(alpha);
                    eta_coef_table.Add(alpha, eta_table);
                }

                BesselYEtaTable eta = eta_table;

                if (!xi_coef_table.TryGetValue(alpha, out BesselYXiTable xi_table)) {
                    xi_table = new BesselYXiTable(alpha, eta);
                    xi_coef_table.Add(alpha, xi_table);
                }

                BesselYXiTable xi = xi_table;

                if (!phi_coef_table.TryGetValue(alpha, out BesselJPhiTable phi_table)) {
                    phi_table = new BesselJPhiTable(alpha);
                    phi_coef_table.Add(alpha, phi_table);
                }

                BesselJPhiTable phi = phi_table;

                Complex f0 = 1e-256, f1 = 0d, lambda = 0d;
                Complex se = 0d, sxo = 0d, sxe = 0d;
                Complex v = 1d / z;

                for (int k = m; k >= 1; k--) {
                    if ((k & 1) == 0) {
                        lambda += f0 * phi[k / 2];

                        se += f0 * eta[k / 2];
                        sxe += f0 * xi[k];
                    }
                    else if (k >= 3) {
                        sxo += f0 * xi[k];
                    }

                    (f0, f1) = (Complex.Ldexp(k + alpha, 1) * v * f0 - f1, f0);
                }

                Complex s = Complex.Pow(Complex.Ldexp(v, 1), alpha), sqs = s * s;

                lambda += f0 * phi[0];
                lambda *= s;

                ddouble rcot = 1d / ddouble.TanPI(alpha), rgamma = ddouble.Gamma(1d + alpha), rsqgamma = rgamma * rgamma;
                Complex r = Complex.Ldexp(ddouble.RcpPI * sqs, 1);
                Complex p = sqs * rsqgamma * ddouble.RcpPI;

                Complex eta0 = ddouble.Abs(alpha) > MillerBwdBesselYEps
                    ? rcot - p / alpha
                    : BesselYEta0Eps(alpha, z);

                Complex xi0 = -Complex.Ldexp(v, 1) * p;
                Complex xi1 = ddouble.Abs(alpha) > MillerBwdBesselYEps
                    ? rcot + p * (alpha * (alpha + 1d) + 1d) / (alpha * (alpha - 1d))
                    : BesselYXi1Eps(alpha, z);

                Complex y0 = r * se + eta0 * f0;
                Complex y1 = r * (3d * alpha * v * sxe + sxo) + xi0 * f0 + xi1 * f1;

                if (n == 0) {
                    Complex yn = y0 / lambda;

                    return yn;
                }
                if (n == 1) {
                    Complex yn = y1 / lambda;

                    return yn;
                }
                if (n >= 0) {
                    for (int k = 1; k < n; k++) {
                        (y1, y0) = (Complex.Ldexp(k + alpha, 1) * v * y1 - y0, y1);
                    }

                    Complex yn = y1 / lambda;

                    return yn;
                }
                else {
                    for (int k = 0; k > n; k--) {
                        (y0, y1) = (Complex.Ldexp(k + alpha, 1) * v * y0 - y1, y0);
                    }

                    Complex yn = y0 / lambda;

                    return yn;
                }
            }

            private static Complex BesselYEta0Eps(ddouble alpha, Complex z) {
                Complex lnz = Complex.Log(z), lnhalfz = Complex.Log(Complex.Ldexp(z, -1));
                ddouble pi = ddouble.PI, sqpi = pi * pi;
                ddouble ln2 = ddouble.Ln2, sqln2 = ln2 * ln2, cbln2 = sqln2 * ln2, qdln2 = sqln2 * sqln2;
                ddouble g = ddouble.EulerGamma;

                Complex r0 = lnhalfz + g;
                Complex r1 =
                    (-sqln2 + lnz * (ln2 * 2d - lnz)) * 4d
                    - sqpi
                    - g * (lnhalfz * 2d + g) * 4d;
                Complex r2 =
                    (-cbln2 + lnz * (sqln2 * 3d + lnz * (ln2 * -3d + lnz))) * 4d
                    + ddouble.Zeta3 * 2d
                    + sqpi * (lnhalfz + g)
                    + g * ((sqln2 + lnz * (ln2 * -2d + lnz)) * 3d + g * (lnhalfz * 3d + g)) * 4d;
                Complex r3 =
                    (-qdln2 + lnz * (cbln2 * 4d + lnz * (sqln2 * -6d + lnz * (ln2 * 4d - lnz)))) * 16d
                    - ddouble.Zeta3 * (lnhalfz + g) * 32d
                    - sqpi * ((sqln2 + lnz * (-ln2 * 2d + lnz) + g * (lnhalfz * 2d + g)) * 8d + sqpi)
                    + g * ((cbln2 + lnz * (sqln2 * -3d + lnz * (ln2 * 3d - lnz))) * 4d
                    + g * ((sqln2 + lnz * (ln2 * -2d + lnz)) * -6d
                    + g * (lnhalfz * -4d
                    - g))) * 16d;

                Complex eta0 = (r0 * 48d + alpha * (r1 * 12d + alpha * (r2 * 8d + alpha * r3))) / (24d * ddouble.PI);

                return eta0;
            }

            static Complex BesselYXi1Eps(ddouble alpha, Complex z) {
                Complex lnz = Complex.Log(z), lnhalfz = Complex.Log(Complex.Ldexp(z, -1)), lnxm1 = lnz - 1, lnhalfxm1 = lnhalfz - 1;
                ddouble pi = ddouble.PI, sqpi = pi * pi;
                ddouble ln2 = ddouble.Ln2, sqln2 = ln2 * ln2, cbln2 = sqln2 * ln2, qdln2 = sqln2 * sqln2;
                ddouble g = ddouble.EulerGamma;

                Complex r0 = lnhalfxm1 + g;
                Complex r1 =
                    (-sqln2 + ln2 * lnxm1 * 2d + lnz * (2 - lnz)) * 4d
                    - sqpi
                    - g * (lnhalfxm1 * 2d + g) * 4d
                    - 6d;
                Complex r2 =
                    -cbln2 * 4d + sqln2 * lnxm1 * 12d + lnz * (18d + lnz * (-12d + lnz * 4d))
                    + ln2 * (lnz * (2d - lnz) * 12d - 18d)
                    + ddouble.Zeta3 * 2d
                    + sqpi * (lnhalfxm1 + g)
                    + g * ((sqln2 - ln2 * lnxm1 * 2d + lnz * (-2d + lnz)) * 12d + 18d
                    + g * (lnhalfxm1 * 12d
                    + g * 4d))
                    - 9d;
                Complex r3 =
                    -qdln2 * 16d
                    + cbln2 * lnxm1 * 64d
                    + sqln2 * (lnz * (2d - lnz) * 96d - 144d)
                    + ln2 * (lnz * (9d + lnz * (-6d + lnz * 2d)) * 32d - 144d)
                    + lnz * (9d + lnz * (-9d + lnz * (4d - lnz))) * 16d
                    + ddouble.Zeta3 * (lnhalfxm1 + g) * -32d
                    + sqpi * ((-sqln2 + ln2 * lnxm1 * 2d + lnz * (2d - lnz) - g * (lnhalfxm1 * 2d + g)) * 8d - 12d - sqpi)
                    + g * ((cbln2 - sqln2 * lnxm1 * 3d) * 64d + ln2 * (lnz * (-2d + lnz) * 192d + 288d) + lnz * (-9d + lnz * (6d - lnz * 2d)) * 32d + 144d
                    + g * ((-sqln2 + ln2 * lnxm1 * 2d + lnz * (2d - lnz)) * 96d - 144d
                    + g * (lnhalfxm1 * -64d
                    - g * 16d)))
                    - 72d;

                Complex xi1 = (r0 * 48d + alpha * (r1 * 12d + alpha * (r2 * 8d + alpha * r3))) / (24d * ddouble.PI);

                return xi1;
            }

            private static Complex BesselY0Kernel(Complex z, int m) {
                Debug.Assert(m >= 2 && (m & 1) == 0);

                if (!eta_coef_table.TryGetValue(0, out BesselYEtaTable eta_table)) {
                    eta_table = new BesselYEtaTable(0);
                    eta_coef_table.Add(0, eta_table);
                }

                BesselYEtaTable eta = eta_table;

                Complex f0 = 1e-256, f1 = 0d, lambda = 0d;
                Complex se = 0d;
                Complex v = 1d / z;

                for (int k = m; k >= 1; k--) {
                    if ((k & 1) == 0) {
                        lambda += f0;

                        se += f0 * eta[k / 2];
                    }

                    (f0, f1) = (2 * k * v * f0 - f1, f0);
                }

                lambda = Complex.Ldexp(lambda, 1) + f0;

                Complex y0 = Complex.Ldexp((se + f0 * (Complex.Log(Complex.Ldexp(z, -1)) + ddouble.EulerGamma)) / (ddouble.PI * lambda), 1);

                return y0;
            }

            private static Complex BesselY1Kernel(Complex z, int m) {
                Debug.Assert(m >= 2 && (m & 1) == 0);

                if (!xi_coef_table.ContainsKey(0)) {
                    if (!eta_coef_table.ContainsKey(0)) {
                        eta_coef_table.Add(0, new BesselYEtaTable(0));
                    }

                    xi_coef_table.Add(0, new BesselYXiTable(0, eta_coef_table[0]));
                }

                BesselYXiTable xi = xi_coef_table[0];

                Complex f0 = 1e-256, f1 = 0d, lambda = 0d;
                Complex sx = 0d;
                Complex v = 1d / z;

                for (int k = m; k >= 1; k--) {
                    if ((k & 1) == 0) {
                        lambda += f0;
                    }
                    else if (k >= 3) {
                        sx += f0 * xi[k];
                    }

                    (f0, f1) = (2 * k * v * f0 - f1, f0);
                }

                lambda = Complex.Ldexp(lambda, 1) + f0;

                Complex y1 = Complex.Ldexp((sx - v * f0 + (Complex.Log(Complex.Ldexp(z, -1)) + ddouble.EulerGamma - 1d) * f1) / (lambda * ddouble.PI), 1);

                return y1;
            }

            private static Complex BesselIKernel(int n, Complex z, int m, bool scale = false) {
                Debug.Assert(m >= 2 && (m & 1) == 0 && n < m);

                n = int.Abs(n);

                if (n == 0) {
                    return BesselI0Kernel(z, m, scale);
                }
                if (n == 1) {
                    return BesselI1Kernel(z, m, scale);
                }

                Complex f0 = 1e-256, f1 = 0d, lambda = 0d, fn = 0d;
                Complex v = 1d / z;

                for (int k = m; k >= 1; k--) {
                    lambda += f0;

                    (f0, f1) = (2 * k * v * f0 + f1, f0);

                    if (k - 1 == n) {
                        fn = f0;
                    }
                }

                lambda = Complex.Ldexp(lambda, 1) + f0;

                Complex yn = fn / lambda;

                if (!scale) {
                    yn *= Complex.Exp(z);
                }

                return yn;
            }

            public static Complex BesselIKernel(ddouble nu, Complex z, int m, bool scale = false) {
                int n = (int)ddouble.Floor(nu);
                ddouble alpha = nu - n;

                if (alpha == 0d) {
                    return BesselIKernel(n, z, m, scale);
                }

                Debug.Assert(m >= 2 && (m & 1) == 0 && n < m);

                if (!psi_coef_table.TryGetValue(alpha, out BesselIPsiTable psi_table)) {
                    psi_table = new BesselIPsiTable(alpha);
                    psi_coef_table.Add(alpha, psi_table);
                }

                BesselIPsiTable psi = psi_table;

                Complex g0 = 1e-256, g1 = 0d, lambda = 0d;
                Complex v = 1d / z;

                if (n >= 0) {
                    Complex gn = 0d;

                    for (int k = m; k >= 1; k--) {
                        lambda += g0 * psi[k];

                        (g0, g1) = (Complex.Ldexp(k + alpha, 1) * v * g0 + g1, g0);

                        if (k - 1 == n) {
                            gn = g0;
                        }
                    }

                    lambda += g0 * psi[0];
                    lambda *= Complex.Pow(Complex.Ldexp(v, 1), alpha);

                    Complex yn = gn / lambda;

                    if (!scale) {
                        yn *= Complex.Exp(z);
                    }

                    return yn;
                }
                else {
                    for (int k = m; k >= 1; k--) {
                        lambda += g0 * psi[k];

                        (g0, g1) = (Complex.Ldexp(k + alpha, 1) * v * g0 + g1, g0);
                    }

                    lambda += g0 * psi[0];
                    lambda *= Complex.Pow(Complex.Ldexp(v, 1), alpha);

                    for (int k = 0; k > n; k--) {
                        (g0, g1) = (Complex.Ldexp(k + alpha, 1) * v * g0 + g1, g0);
                    }

                    Complex yn = g0 / lambda;

                    if (!scale) {
                        yn *= Complex.Exp(z);
                    }

                    return yn;
                }
            }

            private static Complex BesselI0Kernel(Complex z, int m, bool scale = false) {
                Debug.Assert(m >= 2 && (m & 1) == 0);

                Complex g0 = 1e-256, g1 = 0d, lambda = 0d;
                Complex v = 1d / z;

                for (int k = m; k >= 1; k--) {
                    lambda += g0;

                    (g0, g1) = (2 * k * v * g0 + g1, g0);
                }

                lambda = Complex.Ldexp(lambda, 1) + g0;

                Complex y0 = g0 / lambda;

                if (!scale) {
                    y0 *= Complex.Exp(z);
                }

                return y0;
            }

            private static Complex BesselI1Kernel(Complex z, int m, bool scale = false) {
                Debug.Assert(m >= 2 && (m & 1) == 0);

                Complex g0 = 1e-256, g1 = 0d, lambda = 0d;
                Complex v = 1d / z;

                for (int k = m; k >= 1; k--) {
                    lambda += g0;

                    (g0, g1) = (2 * k * v * g0 + g1, g0);
                }

                lambda = Complex.Ldexp(lambda, 1) + g0;

                Complex y1 = g1 / lambda;

                if (!scale) {
                    y1 *= Complex.Exp(z);
                }

                return y1;
            }

            private class BesselJPhiTable {
                private readonly ddouble alpha;
                private readonly List<ddouble> table = [];

                private ddouble g;

                public BesselJPhiTable(ddouble alpha) {
                    Debug.Assert(alpha > 0d && alpha < 1d, nameof(alpha));

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
                    Debug.Assert(alpha > 0d && alpha < 1d, nameof(alpha));

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
                    Debug.Assert(alpha >= 0d && alpha < 1d, nameof(alpha));

                    this.alpha = alpha;
                    this.table.Add(ddouble.NaN);

                    if (alpha > 0d) {
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
                        if (alpha > 0d) {
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
                    Debug.Assert(alpha >= 0d && alpha < 1d, nameof(alpha));

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
                        if (alpha > 0d) {
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
                cds_coef_table.Add(0, Array.AsReadOnly(YoshidaPadeCoefM36.Nu0.Reverse().ToArray()));
                cds_coef_table.Add(1, Array.AsReadOnly(YoshidaPadeCoefM36.Nu1.Reverse().ToArray()));

                List<ReadOnlyCollection<ddouble>> es = [];

                for (int i = 0; i < YoshidaPadeCoefM36.Ess.Length; i++) {
                    es.Add(Array.AsReadOnly(YoshidaPadeCoefM36.Ess[i].ToArray()));
                }

                ess_coef_table = Array.AsReadOnly(es.ToArray());
            }

            public static Complex BesselK(ddouble nu, Complex z, bool scale = false) {
                if (nu < 2d) {
                    if (!cds_coef_table.TryGetValue(nu, out ReadOnlyCollection<(ddouble c, ddouble s)> cds_table)) {
                        cds_table = Table(nu);
                        cds_coef_table.Add(nu, cds_table);
                    }

                    ReadOnlyCollection<(ddouble, ddouble)> cds = cds_table;

                    Complex y = Value(z, cds, scale);

                    return y;
                }
                else {
                    int n = (int)ddouble.Floor(nu);
                    ddouble alpha = nu - n;

                    Complex y0 = BesselK(alpha, z, scale);
                    Complex y1 = BesselK(alpha + 1d, z, scale);

                    Complex v = 1d / z;

                    for (int k = 1; k < n; k++) {
                        (y1, y0) = (Complex.Ldexp(k + alpha, 1) * v * y1 + y0, y1);
                    }

                    return y1;
                }
            }

            private static Complex Value(Complex z, ReadOnlyCollection<(ddouble c, ddouble d)> cds, bool scale = false) {
                Complex t = 1d / z;
                (Complex sc, Complex sd) = cds[0];

                for (int i = 1; i < cds.Count; i++) {
                    (ddouble c, ddouble d) = cds[i];

                    sc = sc * t + c;
                    sd = sd * t + d;
                }

                Complex y = Complex.Sqrt(Complex.Ldexp(t * ddouble.PI, -1)) * sc / sd;

                if (!scale) {
                    y *= Complex.Exp(-z);
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
            public static Complex BesselYPowerSeries(ddouble nu, Complex z) {
                int n = (int)ddouble.Round(nu);
                ddouble alpha = nu - n;

                Complex y0 = PowerSeries.BesselY(n, z);
                Complex y1 = PowerSeries.BesselY(n + ddouble.Sign(alpha) * InterpolationThreshold, z);
                Complex y2 = PowerSeries.BesselY(n + ddouble.Sign(alpha) * InterpolationThreshold * 1.5, z);
                Complex y3 = PowerSeries.BesselY(n + ddouble.Sign(alpha) * InterpolationThreshold * 2, z);

                ddouble t = ddouble.Abs(alpha) / InterpolationThreshold;
                Complex y = Interpolate(t, y0, y1, y2, y3);

                return y;
            }

            public static Complex BesselYMillerBackward(ddouble nu, Complex z) {
                int n = (int)ddouble.Round(nu);
                ddouble alpha = nu - n;

                Complex y0 = MillerBackward.BesselY(n, z);
                Complex y1 = MillerBackward.BesselY(n + ddouble.Sign(alpha) * InterpolationThreshold, z);
                Complex y2 = MillerBackward.BesselY(n + ddouble.Sign(alpha) * InterpolationThreshold * 1.5, z);
                Complex y3 = MillerBackward.BesselY(n + ddouble.Sign(alpha) * InterpolationThreshold * 2, z);

                ddouble t = ddouble.Abs(alpha) / InterpolationThreshold;
                Complex y = Interpolate(t, y0, y1, y2, y3);

                return y;
            }

            public static Complex BesselKPowerSeries(ddouble nu, Complex x) {
                int n = (int)ddouble.Round(nu);
                ddouble alpha = nu - n;

                Complex y0 = PowerSeries.BesselK(n, x);
                Complex y1 = PowerSeries.BesselK(n + ddouble.Sign(alpha) * InterpolationThreshold, x);
                Complex y2 = PowerSeries.BesselK(n + ddouble.Sign(alpha) * InterpolationThreshold * 1.5, x);
                Complex y3 = PowerSeries.BesselK(n + ddouble.Sign(alpha) * InterpolationThreshold * 2, x);

                ddouble t = ddouble.Abs(alpha) / InterpolationThreshold;
                Complex y = Interpolate(t, y0, y1, y2, y3);

                return y;
            }

            private static Complex Interpolate(ddouble t, Complex y0, Complex y1, Complex y2, Complex y3) {
                return y0 + (
                    -(13d + t * (-9d + t * 2d)) / 6d * y0
                    + (6d + t * (-7d + t * 2d)) * y1
                    - (16d + t * (-24d + t * 8d)) / 3d * y2
                    + (3d + t * (-5d + t * 2d)) / 2d * y3) * t;
            }
        }
    }
}