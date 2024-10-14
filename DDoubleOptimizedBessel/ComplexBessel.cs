using DoubleDouble;
using DoubleDoubleComplex;
using System.Collections.ObjectModel;
using System.Diagnostics;
using static DDoubleOptimizedBessel.ComplexBesselUtil;

namespace DDoubleOptimizedBessel {

    public static class ComplexBessel {

        public static Complex BesselJ(ddouble nu, Complex z) {
            CheckNu(nu);

            if (!Complex.IsFinite(z)) {
                return Complex.NaN;
            }

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
                return (SinCosPICache.CosPI(nu / 2d), SinCosPICache.SinPI(nu / 2d)) * MillerBackward.BesselI(nu, (z.I, z.R)).Conj;
            }
        }

        public static Complex BesselY(ddouble nu, Complex z) {
            CheckNu(nu);

            if (!Complex.IsFinite(z)) {
                return Complex.NaN;
            }

            if (ddouble.IsNegative(z.I)) {
                return BesselY(nu, z.Conj).Conj;
            }

            if (ddouble.IsNegative(z.R)) {
                return (SinCosPICache.CosPI(nu), -SinCosPICache.SinPI(nu)) * BesselY(nu, -z)
                        + (0d, 2d * SinCosPICache.CosPI(nu)) * BesselJ(nu, -z);
            }

            if (z.Magnitude >= HankelThreshold) {
                return Limit.BesselY(nu, z);
            }
            else if (z.R <= PowerSeriesThreshold(nu, z.I) - BesselJYPowerseriesBias) {
                if (NearlyInteger(nu, out int n) || ddouble.Abs(n - nu) >= BesselYForcedMillerBackwardThreshold) {
                    return PowerSeries.BesselY(nu, z);
                }
                else if (z.I <= MillerBackwardThreshold / 2d) {
                    return MillerBackward.BesselY(nu, z);
                }
                else if (ddouble.Abs(n - nu) >= BesselYInterpolationDelta) {
                    return PowerSeries.BesselY(nu, z);
                }
                else {
                    return Interpolation.BesselY(nu, z);
                }
            }
            else if (z.I <= MillerBackwardThreshold) {
                return MillerBackward.BesselY(nu, z);
            }
            else {
                Complex c = (SinCosPICache.CosPI(nu / 2d), SinCosPICache.SinPI(nu / 2d));

                Complex bi = (z.R <= PowerSeriesThreshold(nu, z.I))
                    ? PowerSeries.BesselI(nu, (z.I, z.R))
                    : MillerBackward.BesselI(nu, (z.I, z.R));

                Complex bk = YoshidaPade.BesselK(nu, (z.I, z.R));

                Complex y = Complex.ImaginaryOne * c * bi.Conj - 2d * ddouble.RcpPI * (c * bk).Conj;

                return y;
            }
        }

        public static Complex BesselI(ddouble nu, Complex z) {
            CheckNu(nu);

            if (!Complex.IsFinite(z)) {
                return Complex.NaN;
            }

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

            if (!Complex.IsFinite(z)) {
                return Complex.NaN;
            }

            nu = ddouble.Abs(nu);

            if (ddouble.IsNegative(z.I)) {
                return BesselK(nu, z.Conj).Conj;
            }

            if (ddouble.IsNegative(z.R)) {
                return (SinCosPICache.CosPI(nu), -SinCosPICache.SinPI(nu)) * BesselK(nu, -z)
                        - (0d, ddouble.PI) * BesselI(nu, -z);
            }

            if (z.Magnitude >= HankelThreshold) {
                return Limit.BesselK(nu, z);
            }
            else if (z.Magnitude <= BesselKNearZeroThreshold) {
                if (NearlyInteger(nu, out int n) || ddouble.Abs(n - nu) >= BesselKInterpolationDelta) {
                    return PowerSeries.BesselK(nu, z);
                }
                else {
                    return Interpolation.BesselK(nu, z);
                }
            }
            else if (z.R >= BesselKPadeThreshold) {
                return YoshidaPade.BesselK(nu, z);
            }
            else {
                Complex c = (SinCosPICache.CosPI(nu / 2d), -SinCosPICache.SinPI(nu / 2d));

                Complex bi = (z.I <= PowerSeriesThreshold(nu, z.R))
                    ? PowerSeries.BesselI(nu, z)
                    : MillerBackward.BesselI(nu, z);

                Complex by = (
                    (z.I <= PowerSeriesThreshold(nu, z.R) - BesselJYPowerseriesBias)
                    ? (
                        (NearlyInteger(nu, out int n) || ddouble.Abs(n - nu) >= BesselYForcedMillerBackwardThreshold)
                         ? PowerSeries.BesselY(nu, (z.I, z.R))
                         : (z.R <= MillerBackwardThreshold / 2d)
                         ? MillerBackward.BesselY(nu, (z.I, z.R))
                         : (ddouble.Abs(n - nu) >= BesselYInterpolationDelta)
                         ? PowerSeries.BesselY(nu, (z.I, z.R))
                         : Interpolation.BesselY(nu, (z.I, z.R))
                    )
                    : MillerBackward.BesselY(nu, (z.I, z.R))
                );

                Complex y = c * (-Complex.ImaginaryOne * c * bi - by.Conj) * ddouble.PI / 2d;

                return y;
            }
        }

        public static Complex HankelH1(ddouble nu, Complex z) {
            Complex y = BesselJ(nu, z) + Complex.ImaginaryOne * BesselY(nu, z);

            return y;
        }

        public static Complex HankelH2(ddouble nu, Complex z) {
            Complex y = BesselJ(nu, z) - Complex.ImaginaryOne * BesselY(nu, z);

            return y;
        }
    }

    static class ComplexBesselUtil {
        public const int MaxN = 16;
        public static readonly double Eps = double.ScaleB(1, -105);
        public static readonly int NearZeroExponent = -950;
        public static readonly double BesselYForcedMillerBackwardThreshold = double.ScaleB(1, -8);
        public static readonly double BesselYInterpolationDelta = double.ScaleB(1, -16);
        public static readonly double BesselKInterpolationDelta = double.ScaleB(1, -8);
        public const double HankelThreshold = 38.875, MillerBackwardThreshold = 6;
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
                    Complex y = BesselJKernel(nu, z, terms: 42);

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

            public static Complex BesselI(ddouble nu, Complex z) {
                Debug.Assert(ddouble.IsPositive(z.R));
                Debug.Assert(ddouble.IsPositive(z.I));

                if (ddouble.IsNegative(nu) && NearlyInteger(nu, out _)) {
                    Complex y = BesselI(-nu, z);

                    return y;
                }
                else {
                    Complex y = BesselIKernel(nu, z, terms: 42);

                    return y;
                }
            }

            public static Complex BesselK(ddouble nu, Complex z) {
                Debug.Assert(ddouble.IsPositive(z.R));
                Debug.Assert(ddouble.IsPositive(z.I));

                if (NearlyInteger(nu, out int n)) {
                    Complex y = BesselKKernel(n, z, terms: 27);

                    return y;
                }
                else {
                    Complex y = BesselKKernel(nu, z, terms: 30);

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

                if (!Complex.IsFinite(u) || Complex.IsZero(z2)) {
                    if (ddouble.IsZero(nu)) {
                        return 1d;
                    }
                    if (NearlyInteger(nu, out _) || ddouble.IsPositive(nu)) {
                        return 0d;
                    }

                    return Complex.NaN;
                }

                for (int k = 0; k <= terms; k++) {
                    c = SeriesUtil.Add(c, u * r[k], 1d, -z2 * d[k], out bool convergence);

                    if (convergence) {
                        break;
                    }

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
                Complex p = Complex.IsZero(cos) ? 0d : Complex.Pow(z, Complex.Ldexp(nu, 1)) * cos;
                Complex s = Complex.Ldexp(Complex.Pow(Complex.Ldexp(z, 1), nu), 2);

                Complex z2 = z * z, z4 = z2 * z2;

                if (!Complex.IsFinite(p) || !Complex.IsFinite(s) || Complex.ILogB(s) < NearZeroExponent || Complex.IsZero(z2)) {
                    return Complex.NaN;
                }

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

                if (!Complex.IsFinite(h) || Complex.IsZero(z2)) {
                    return Complex.NaN;
                }

                Complex c = 0d, u = Complex.Ldexp(ddouble.RcpPI, 1);

                for (int k = 0; k <= terms; k++) {
                    Complex s = u * r[k], t = h - ddouble.HarmonicNumber(2 * k);
                    c = SeriesUtil.Add(c, s * t, 1d, -z2 * d[k], out bool convergence1);
                    c = SeriesUtil.Add(c, s, z2 * q[k], out bool convergence2);

                    if (convergence1 && convergence2 && Complex.ILogB(t) >= -4) {
                        break;
                    }

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

                if (!Complex.IsFinite(h) || Complex.IsZero(z2)) {
                    return Complex.NaN;
                }

                Complex c = -2d / (z * ddouble.PI), u = z / Complex.Ldexp(ddouble.PI, 1);

                for (int k = 0; k <= terms; k++) {
                    Complex s = u * r[k], t = h - ddouble.HarmonicNumber(2 * k) - ddouble.HarmonicNumber(2 * k + 1);
                    c = SeriesUtil.Add(c, s * t, 1d, -z2 * d[k], out bool convergence1);
                    c = SeriesUtil.Add(c, s, z2 * q[k], out bool convergence2);

                    if (convergence1 && convergence2 && Complex.ILogB(t) >= -4) {
                        break;
                    }

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

                Complex h = Complex.Ldexp(Complex.Log(Complex.Ldexp(z, -1)) + ddouble.EulerGamma, 1);

                if (!Complex.IsFinite(h)) {
                    return Complex.NaN;
                }

                Complex c = 0d;
                Complex z2 = z * z, z4 = z2 * z2;
                Complex u = 1d, v = 1d, w = Complex.Ldexp(z2, -2);

                for (int k = 0; k < n; k++) {
                    c += v * f[k];
                    v *= w;
                }
                c /= -v;

                if (!Complex.IsFinite(c) || Complex.IsZero(z2)) {
                    return Complex.NaN;
                }

                for (int k = 0; k <= terms; k++) {
                    Complex s = u * r[k], t = (h - ddouble.HarmonicNumber(2 * k) - ddouble.HarmonicNumber(2 * k + n));
                    c = SeriesUtil.Add(c, s * t, 1d, -z2 * d[k], out bool convergence1);
                    c = SeriesUtil.Add(c, s, z2 * q[k], out bool convergence2);

                    if (convergence1 && convergence2 && Complex.ILogB(t) >= -4) {
                        break;
                    }

                    u *= z4;
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

                if (!Complex.IsFinite(u) || Complex.IsZero(z2)) {
                    if (ddouble.IsZero(nu)) {
                        return 1d;
                    }
                    if (NearlyInteger(nu, out _) || ddouble.IsPositive(nu)) {
                        return 0d;
                    }

                    return Complex.NaN;
                }

                for (int k = 0; k <= terms; k++) {
                    c = SeriesUtil.Add(c, u * r[k], 1d, z2 * d[k], out bool convergence);

                    if (convergence) {
                        break;
                    }

                    u *= z4;

                    if (!Complex.IsFinite(c)) {
                        break;
                    }
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

                if (Complex.ILogB(tp) < NearZeroExponent || Complex.IsZero(z2)) {
                    return Complex.NaN;
                }

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

                if (!Complex.IsFinite(h)) {
                    return Complex.NaN;
                }

                Complex z2 = z * z;

                if (Complex.IsZero(z2)) {
                    return Complex.NaN;
                }

                Complex c = 0d, u = 1d;

                for (int k = 0; k <= terms; k++) {
                    c = SeriesUtil.Add(c, u * r[k], h, ddouble.HarmonicNumber(k), out bool convergence);

                    if (convergence) {
                        break;
                    }

                    u *= z2;
                }

                return c;
            }

            private static Complex BesselK1Kernel(Complex z, int terms) {
                K1CoefTable r = k1_coef_table;
                Complex h = Complex.Log(Complex.Ldexp(z, -1)) + ddouble.EulerGamma;

                if (!Complex.IsFinite(h)) {
                    return Complex.NaN;
                }

                Complex z2 = z * z;

                if (Complex.IsZero(z2)) {
                    return Complex.NaN;
                }

                Complex c = 1d / z, u = Complex.Ldexp(z, -1);

                for (int k = 0; k <= terms; k++) {
                    c = SeriesUtil.Add(c, u * r[k], h, -ddouble.Ldexp(ddouble.HarmonicNumber(k) + ddouble.HarmonicNumber(k + 1), -1), out bool convergence);

                    if (convergence) {
                        break;
                    }

                    u *= z2;
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

        public static class Limit {
            static readonly Dictionary<ddouble, HankelExpansion> table = [];

            public static Complex BesselJ(ddouble nu, Complex z) {
                Debug.Assert(ddouble.IsPositive(z.R));
                Debug.Assert(ddouble.IsPositive(z.I));

                if (!table.TryGetValue(nu, out HankelExpansion hankel)) {
                    hankel = new HankelExpansion(nu);
                    table.Add(nu, hankel);
                }

                if (!Complex.IsFinite(z)) {
                    return Complex.NaN;
                }

                (Complex c_even, Complex c_odd) = hankel.BesselJYCoef(z);

                Complex omega = hankel.Omega(z);

                Complex cos = Complex.Cos(omega), sin = Complex.Sin(omega);

                Complex y = Complex.Sqrt(2d / (ddouble.PI * z)) * (cos * c_even - sin * c_odd);

                return y;
            }

            public static Complex BesselY(ddouble nu, Complex z) {
                Debug.Assert(ddouble.IsPositive(z.R));
                Debug.Assert(ddouble.IsPositive(z.I));

                if (!table.TryGetValue(nu, out HankelExpansion hankel)) {
                    hankel = new HankelExpansion(nu);
                    table.Add(nu, hankel);
                }

                if (!Complex.IsFinite(z)) {
                    return Complex.NaN;
                }

                (Complex c_even, Complex c_odd) = hankel.BesselJYCoef(z);

                Complex omega = hankel.Omega(z);

                Complex cos = Complex.Cos(omega), sin = Complex.Sin(omega);

                Complex y = Complex.Sqrt(2d / (ddouble.PI * z)) * (sin * c_even + cos * c_odd);

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

                Complex y = Complex.Sqrt(1d / (2d * ddouble.PI * z)) * (
                    Complex.Exp(z) * ci -
                    (SinCosPICache.SinPI(nu), -SinCosPICache.CosPI(nu)) * Complex.Exp(-z) * ck
                );

                if (!Complex.IsFinite(z) || Complex.IsZero(y)) {
                    return Complex.NaN;
                }

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

                Complex y = Complex.Sqrt(ddouble.PI / (2d * z)) * Complex.Exp(-z) * c;

                if (!Complex.IsFinite(z) || Complex.IsZero(y)) {
                    return Complex.NaN;
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

                public Complex Omega(Complex z) {
                    Complex omega = z - ddouble.Ldexp(2d * Nu + 1d, -2) * ddouble.PI;

                    return omega;
                }

                public (Complex c_even, Complex c_odd) BesselJYCoef(Complex z, int terms = 35) {
                    Complex v = 1d / (z * z), w = -v;

                    Complex c_even = ACoef(0), c_odd = ACoef(1);

                    for (int k = 1; k <= terms; k++) {
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

                public Complex BesselICoef(Complex z, int terms = 75) {
                    Complex v = 1d / z, w = -v;

                    Complex c = ACoef(0);

                    for (int k = 1; k <= terms; k++) {
                        Complex dc = w * ACoef(k);

                        c += dc;

                        if ((long)Complex.ILogB(c) - Complex.ILogB(dc) >= 106L || Complex.IsZero(dc)) {
                            break;
                        }

                        w *= -v;
                    }

                    return c;
                }

                public Complex BesselKCoef(Complex z, int terms = 58) {
                    Complex v = 1d / z, w = v;

                    Complex c = ACoef(0);

                    for (int k = 1; k <= terms; k++) {
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
            public static readonly int BesselYEpsExponent = -12;

            private static readonly Dictionary<ddouble, BesselJPhiTable> phi_coef_table = [];
            private static readonly Dictionary<ddouble, BesselIPsiTable> psi_coef_table = [];
            private static readonly Dictionary<ddouble, BesselYEtaTable> eta_coef_table = [];
            private static readonly Dictionary<ddouble, BesselYXiTable> xi_coef_table = [];
            private static readonly ReadOnlyCollection<ddouble> eta0_coef = new(MillerBackwardCoef.Eta0.Reverse().ToArray());
            private static readonly ReadOnlyCollection<ddouble> xi1_coef = new(MillerBackwardCoef.Xi1.Reverse().ToArray());

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

            public static Complex BesselI(int n, Complex z) {
                Debug.Assert(ddouble.IsPositive(z.R));
                Debug.Assert(ddouble.IsPositive(z.I));

                int m = BesselIIterM((double)z.R, (double)z.I);

                Complex y = BesselIKernel(n, z, m);

                return y;
            }

            public static Complex BesselI(ddouble nu, Complex z) {
                Debug.Assert(ddouble.IsPositive(z.R));
                Debug.Assert(ddouble.IsPositive(z.I));

                int m = BesselIIterM((double)z.R, (double)z.I);

                if (NearlyInteger(nu, out int n)) {
                    Complex y = BesselIKernel(n, z, m);

                    return y;
                }
                else {
                    Complex y = BesselIKernel(nu, z, m);

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
                else if (n == 0) {
                    return BesselJ0Kernel(z, m);
                }
                else if (n == 1) {
                    return BesselJ1Kernel(z, m);
                }
                else {
                    return BesselJNKernel(n, z, m);
                }
            }

            private static Complex BesselJKernel(ddouble nu, Complex z, int m) {
                int n = (int)ddouble.Floor(nu);
                ddouble alpha = nu - n;

                Debug.Assert(m >= 2 && (m & 1) == 0 && n < m);

                if (!phi_coef_table.TryGetValue(alpha, out BesselJPhiTable phi)) {
                    phi = new BesselJPhiTable(alpha);
                    phi_coef_table.Add(alpha, phi);
                }

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

            private static Complex BesselJNKernel(int n, Complex z, int m) {
                Debug.Assert(m >= 2 && (m & 1) == 0);

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

            private static Complex BesselYKernel(int n, Complex z, int m) {
                Debug.Assert(m >= 2 && (m & 1) == 0 && n < m);

                if (n < 0) {
                    return (n & 1) == 0 ? BesselYKernel(-n, z, m) : -BesselYKernel(-n, z, m);
                }
                else if (n == 0) {
                    return BesselY0Kernel(z, m);
                }
                else if (n == 1) {
                    return BesselY1Kernel(z, m);
                }
                else {
                    return BesselYNKernel(n, z, m);
                }
            }

            private static Complex BesselYKernel(ddouble nu, Complex z, int m) {
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

                Complex xi0 = -Complex.Ldexp(v, 1) * p;

                (Complex eta0, Complex xi1) = ddouble.ILogB(alpha) >= BesselYEpsExponent
                    ? (rcot - p / alpha, rcot + p * (alpha * (alpha + 1d) + 1d) / (alpha * (alpha - 1d)))
                    : BesselYEta0Xi1Eps(alpha, z);

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

            private static (Complex eta0, Complex xi1) BesselYEta0Xi1Eps(ddouble alpha, Complex z) {
                const int N = 7;

                Complex lnz = Complex.Log(z);

                Complex eta0 = 0d, xi1 = 0d;
                for (int i = N, k = 0; i >= 0; i--) {
                    Complex s = eta0_coef[k], t = xi1_coef[k];
                    k++;

                    for (int j = i; j >= 0; j--, k++) {
                        s = eta0_coef[k] + lnz * s;
                        t = xi1_coef[k] + lnz * t;
                    }

                    eta0 = s + alpha * eta0;
                    xi1 = t + alpha * xi1;
                }

                return (eta0, xi1);
            }

            private static Complex BesselY0Kernel(Complex z, int m) {
                Debug.Assert(m >= 2 && (m & 1) == 0);

                if (!eta_coef_table.TryGetValue(0, out BesselYEtaTable eta)) {
                    eta = new BesselYEtaTable(0);
                    eta_coef_table.Add(0, eta);
                }

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

                if (!xi_coef_table.TryGetValue(0, out BesselYXiTable xi)) {
                    if (!eta_coef_table.TryGetValue(0, out BesselYEtaTable eta)) {
                        eta = new BesselYEtaTable(0);
                        eta_coef_table.Add(0, eta);
                    }

                    xi = new BesselYXiTable(0, eta);
                    xi_coef_table.Add(0, xi);
                }

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

            private static Complex BesselYNKernel(int n, Complex z, int m) {
                Debug.Assert(m >= 2 && (m & 1) == 0);

                if (!eta_coef_table.TryGetValue(0, out BesselYEtaTable eta)) {
                    eta = new BesselYEtaTable(0);
                    eta_coef_table.Add(0, eta);
                }
                if (!xi_coef_table.TryGetValue(0, out BesselYXiTable xi)) {
                    xi = new BesselYXiTable(0, eta);
                    xi_coef_table.Add(0, xi);
                }

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

            private static Complex BesselIKernel(int n, Complex z, int m) {
                Debug.Assert(m >= 2 && (m & 1) == 0 && n < m);

                n = int.Abs(n);

                if (n == 0) {
                    return BesselI0Kernel(z, m);
                }
                else if (n == 1) {
                    return BesselI1Kernel(z, m);
                }
                else {
                    return BesselINKernel(n, z, m);
                }
            }

            private static Complex BesselIKernel(ddouble nu, Complex z, int m) {
                int n = (int)ddouble.Floor(nu);
                ddouble alpha = nu - n;

                Debug.Assert(m >= 2 && (m & 1) == 0 && n < m);

                if (!psi_coef_table.TryGetValue(alpha, out BesselIPsiTable psi)) {
                    psi = new BesselIPsiTable(alpha);
                    psi_coef_table.Add(alpha, psi);
                }

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

                    yn *= Complex.Exp(z);

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

                    yn *= Complex.Exp(z);

                    return yn;
                }
            }

            private static Complex BesselI0Kernel(Complex z, int m) {
                Debug.Assert(m >= 2 && (m & 1) == 0);

                Complex g0 = 1e-256, g1 = 0d, lambda = 0d;
                Complex v = 1d / z;

                for (int k = m; k >= 1; k--) {
                    lambda += g0;

                    (g0, g1) = (2 * k * v * g0 + g1, g0);
                }

                lambda = Complex.Ldexp(lambda, 1) + g0;

                Complex y0 = g0 / lambda;

                y0 *= Complex.Exp(z);

                return y0;
            }

            private static Complex BesselI1Kernel(Complex z, int m) {
                Debug.Assert(m >= 2 && (m & 1) == 0);

                Complex g0 = 1e-256, g1 = 0d, lambda = 0d;
                Complex v = 1d / z;

                for (int k = m; k >= 1; k--) {
                    lambda += g0;

                    (g0, g1) = (2 * k * v * g0 + g1, g0);
                }

                lambda = Complex.Ldexp(lambda, 1) + g0;

                Complex y1 = g1 / lambda;

                y1 *= Complex.Exp(z);

                return y1;
            }

            private static Complex BesselINKernel(int n, Complex z, int m) {
                Debug.Assert(m >= 2 && (m & 1) == 0);

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

                yn *= Complex.Exp(z);

                return yn;
            }

            private class BesselJPhiTable {
                private readonly ddouble alpha;
                private readonly List<ddouble> table = [];

                private ddouble g;

                public BesselJPhiTable(ddouble alpha) {
                    Debug.Assert(alpha > -1d && alpha < 1d, nameof(alpha));

                    this.alpha = alpha;

                    ddouble phi0 = ddouble.Gamma(1d + alpha);
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
            }

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
            }

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
            }

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
                cds_coef_table.Add(0, Array.AsReadOnly(YoshidaPadeCoefM36.Nu0.Reverse().ToArray()));
                cds_coef_table.Add(1, Array.AsReadOnly(YoshidaPadeCoefM36.Nu1.Reverse().ToArray()));

                List<ReadOnlyCollection<ddouble>> es = [];

                for (int i = 0; i < YoshidaPadeCoefM36.Ess.Length; i++) {
                    es.Add(Array.AsReadOnly(YoshidaPadeCoefM36.Ess[i].ToArray()));
                }

                ess_coef_table = Array.AsReadOnly(es.ToArray());
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

                Complex y = Complex.Sqrt(Complex.Ldexp(t * ddouble.PI, -1)) * sc / sd;

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

        public static class Interpolation {
            public static Complex BesselY(ddouble nu, Complex z) {
                int n = (int)ddouble.Round(nu);
                ddouble alpha = nu - n;

                Debug.Assert(ddouble.Abs(alpha) <= BesselYInterpolationDelta);

                Complex y0 = PowerSeries.BesselY(n, z);

                if (!Complex.IsFinite(y0)) {
                    return y0;
                }

                ddouble dnu = ddouble.Sign(alpha) * BesselYInterpolationDelta;

                Complex y1 = PowerSeries.BesselY(n + dnu, z);
                Complex y2 = PowerSeries.BesselY(n + dnu * 1.25d, z);
                Complex y3 = PowerSeries.BesselY(n + dnu * 1.5d, z);
                Complex y4 = PowerSeries.BesselY(n + dnu * 1.75d, z);
                Complex y5 = PowerSeries.BesselY(n + dnu * 2d, z);

                if (!Complex.IsFinite(y1) || !Complex.IsFinite(y2) || !Complex.IsFinite(y3) || !Complex.IsFinite(y4) || !Complex.IsFinite(y5)) {
                    return y1;
                }

                ddouble t = ddouble.Abs(alpha) / BesselYInterpolationDelta;
                Complex y = Interpolate(t, y0, y1, y2, y3, y4, y5);

                return y;
            }

            public static Complex BesselK(ddouble nu, Complex z) {
                int n = (int)ddouble.Round(nu);
                ddouble alpha = nu - n;

                Debug.Assert(n >= 0);
                Debug.Assert(ddouble.Abs(alpha) <= BesselKInterpolationDelta);

                Complex y0 = PowerSeries.BesselK(0, z);

                if (!Complex.IsFinite(y0)) {
                    return y0;
                }
                ddouble dnu = BesselKInterpolationDelta;

                Complex y1 = PowerSeries.BesselK(dnu, z);
                Complex y2 = PowerSeries.BesselK(dnu * 1.25d, z);
                Complex y3 = PowerSeries.BesselK(dnu * 1.5d, z);
                Complex y4 = PowerSeries.BesselK(dnu * 1.75d, z);
                Complex y5 = PowerSeries.BesselK(dnu * 2d, z);

                if (!Complex.IsFinite(y1) || !Complex.IsFinite(y2) || !Complex.IsFinite(y3) || !Complex.IsFinite(y4) || !Complex.IsFinite(y5)) {
                    return y1;
                }

                ddouble t = ddouble.Abs(alpha) / BesselKInterpolationDelta;
                Complex k0 = InterpolateEvenConvex(t, y0, y1, y2, y3, y4, y5);

                if (n == 0) {
                    return k0;
                }

                Complex i0 = PowerSeries.BesselI(alpha, z), i1 = PowerSeries.BesselI(alpha + 1d, z);

                Complex k1 = (1d - i1 * k0 * z) / (i0 * z);

                if (n == 1) {
                    return k1;
                }

                Complex v = 1d / z;

                for (int k = 1; k < n; k++) {
                    (k1, k0) = (ddouble.Ldexp(k + alpha, 1) * v * k1 + k0, k1);
                }

                return k1;
            }

            private static Complex Interpolate(ddouble t, Complex y0, Complex y1, Complex y2, Complex y3, Complex y4, Complex y5) {
                Complex y = y0
                    + t * (-743d * y0 + 14700d * y1 - 37632d * y2 + 39200d * y3 - 19200d * y4 + 3675d * y5
                    + t * (1035d * y0 - 37310d * y1 + 103040d * y2 - 112560d * y3 + 56960d * y4 - 11165d * y5
                    + t * (-710d * y0 + 35140d * y1 - 103040d * y2 + 118160d * y3 - 62080d * y4 + 12530d * y5
                    + t * (240d * y0 - 14560d * y1 + 44800d * y2 - 53760d * y3 + 29440d * y4 - 6160d * y5
                    + t * (-32d * y0 + 2240d * y1 - 7168d * y2 + 8960d * y3 - 5120d * y4 + 1120d * y5))))) / 210d;

                return y;
            }

            private static Complex InterpolateEvenConvex(ddouble t, Complex y0, Complex y1, Complex y2, Complex y3, Complex y4, Complex y5) {
                ddouble t2 = t * t;

                Complex y = y0
                    + t2 * (-151028163d * y0 + 561834000d * y1 - 708083712d * y2 + 395136000d * y3 - 110592000d * y4 + 12733875d * y5
                    + t2 * (150533955d * y0 - 933192260d * y1 + 1431019520d * y2 - 875831040d * y3 + 258170880d * y4 - 30701055d * y5
                    + t2 * (-70594524d * y0 + 556941840d * y1 - 962174976d * y2 + 658748160d * y3 - 209018880d * y4 + 26098380d * y5
                    + t2 * (15649920d * y0 - 141872640d * y1 + 264929280d * y2 - 198696960d * y3 + 69304320d * y4 - 9313920d * y5
                    + t2 * (-1317888d * y0 + 13045760d * y1 - 25690112d * y2 + 20643840d * y3 - 7864320d * y4 + 1182720d * y5))))) / 56756700d;

                return y;
            }
        }
    }
}