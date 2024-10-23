using MultiPrecision;
using MultiPrecisionComplex;
using System.Diagnostics;

namespace ComplexBessel {
    public static class Recurrence {
        public static Complex<Pow2.N4> BesselJ(MultiPrecision<Pow2.N4> nu, Complex<Pow2.N4> z) {
            Debug.Assert(nu >= BesselUtil<Pow2.N4>.DirectMaxN || nu <= -BesselUtil<Pow2.N4>.DirectMaxN);

            MultiPrecision<Pow2.N4> nu_abs = MultiPrecision<Pow2.N4>.Abs(nu);
            int n = (int)MultiPrecision<Pow2.N4>.Floor(nu_abs);
            MultiPrecision<Pow2.N4> alpha = nu_abs - n;

            Complex<Pow2.N4> v = 1d / z;

            if (MultiPrecision<Pow2.N4>.IsPositive(nu)) {
                (Complex<Pow2.N4> a0, Complex<Pow2.N4> b0, Complex<Pow2.N4> a1, Complex<Pow2.N4> b1) = (1d, 0d, 0d, 1d);

                Complex<Pow2.N4> r = Complex<Pow2.N4>.Ldexp(nu_abs * v, 1);
                (a0, b0, a1, b1) = (a1, b1, r * a1 + a0, r * b1 + b0);

                Complex<Pow2.N4> s = 1d;

                for (int i = 1; i <= 1024; i++) {
                    r = Complex<Pow2.N4>.Ldexp((nu_abs + i) * v, 1);

                    (a0, b0, a1, b1) = (a1, b1, r * a1 - a0, r * b1 - b0);
                    s = a1 / b1;

                    long exp = long.Max(a1.Exponent, b1.Exponent);
                    if (exp <= int.MinValue) {
                        return Complex<Pow2.N4>.NaN;
                    }

                    (a0, a1, b0, b1) = (Complex<Pow2.N4>.Ldexp(a0, -exp), Complex<Pow2.N4>.Ldexp(a1, -exp), Complex<Pow2.N4>.Ldexp(b0, -exp), Complex<Pow2.N4>.Ldexp(b1, -exp));

                    if (i > 0 && (i & 3) == 0) {
                        Complex<Pow2.N4> r0 = a0 * b1, r1 = a1 * b0;
                        if (!((r0 - r1).Magnitude > MultiPrecision<Pow2.N4>.Min(r0.Magnitude, r1.Magnitude) * 1e-36)) {
                            break;
                        }
                    }
                }

                long exp_sum = 0;
                bool s_overone = s.Magnitude > 1d;
                (Complex<Pow2.N4> j0, Complex<Pow2.N4> j1) = s_overone ? (Complex<Pow2.N4>.One, 1d / s) : (s, 1d);

                for (int k = n - 1; k >= BesselUtil<Pow2.N4>.DirectMaxN - 1; k--) {
                    (j1, j0) = (MultiPrecision<Pow2.N4>.Ldexp(k + alpha, 1) * v * j1 - j0, j1);

                    long j0_exp = j0.Exponent, j1_exp = j1.Exponent;
                    if (long.Sign(j0_exp) * long.Sign(j1_exp) > 0) {
                        long exp = j0_exp > 0 ? long.Max(j0_exp, j1_exp) : long.Min(j0_exp, j1_exp);
                        exp_sum += exp;
                        (j0, j1) = (Complex<Pow2.N4>.Ldexp(j0, -exp), Complex<Pow2.N4>.Ldexp(j1, -exp));
                    }
                }

                Complex<Pow2.N4> y = Complex<Pow2.N4>.Ldexp(
                    j0.Magnitude >= j1.Magnitude
                    ? BesselN4.BesselJ(alpha + (BesselUtil<Pow2.N4>.DirectMaxN - 1), z) / j0
                    : BesselN4.BesselJ(alpha + (BesselUtil<Pow2.N4>.DirectMaxN - 2), z) / j1,
                    long.Clamp(-exp_sum, int.MinValue, int.MaxValue)
                ) * (s_overone ? 1d : s);

                return y;
            }
            else {
                if (BesselUtil<Pow2.N4>.NearlyInteger(nu, out int near_n)) {
                    return (near_n & 1) == 0 ? BesselJ(-near_n, z) : -BesselJ(-near_n, z);
                }

                return (SinCosPICache<Pow2.N4>.CosPI(nu / 2), SinCosPICache<Pow2.N4>.SinPI(nu / 2)) * BesselI(nu, (z.I, z.R)).Conj;
            }
        }

        public static Complex<Pow2.N4> BesselY(MultiPrecision<Pow2.N4> nu, Complex<Pow2.N4> z) {
            Debug.Assert(nu >= BesselUtil<Pow2.N4>.DirectMaxN || nu <= -BesselUtil<Pow2.N4>.DirectMaxN);

            if (BesselUtil<Pow2.N4>.NearlyInteger(nu + 0.5d, out int near_n)) {
                return (near_n & 1) == 0 ? BesselJ(-nu, z) : -BesselJ(-nu, z);
            }

            MultiPrecision<Pow2.N4> nu_abs = MultiPrecision<Pow2.N4>.Abs(nu);

            Complex<Pow2.N4> c = (SinCosPICache<Pow2.N4>.CosPI(nu_abs / 2), SinCosPICache<Pow2.N4>.SinPI(nu_abs / 2));
            Complex<Pow2.N4> bi = BesselI(nu_abs, (z.I, z.R));
            Complex<Pow2.N4> bk = BesselK(nu_abs, (z.I, z.R));

            if (MultiPrecision<Pow2.N4>.IsPositive(nu)) {
                Complex<Pow2.N4> y = (0, 1) * c * bi.Conj - 2 * MultiPrecision<Pow2.N4>.RcpPI * (c * bk).Conj;

                return y;
            }
            else {
                Complex<Pow2.N4> y = (0, 1) * (c * bi).Conj
                    + 2 * MultiPrecision<Pow2.N4>.RcpPI * ((0, 1) * c.Conj * SinCosPICache<Pow2.N4>.SinPI(nu_abs) - c) * bk.Conj;

                return y;
            }
        }

        public static Complex<Pow2.N4> BesselI(MultiPrecision<Pow2.N4> nu, Complex<Pow2.N4> z) {
            Debug.Assert(nu >= BesselUtil<Pow2.N4>.DirectMaxN || nu <= -BesselUtil<Pow2.N4>.DirectMaxN);

            MultiPrecision<Pow2.N4> nu_abs = MultiPrecision<Pow2.N4>.Abs(nu);
            int n = (int)MultiPrecision<Pow2.N4>.Floor(nu_abs);
            MultiPrecision<Pow2.N4> alpha = nu_abs - n;

            Complex<Pow2.N4> v = 1d / z;

            (Complex<Pow2.N4> a0, Complex<Pow2.N4> b0, Complex<Pow2.N4> a1, Complex<Pow2.N4> b1) = (1d, 0d, 0d, 1d);
            Complex<Pow2.N4> s = 1d;

            for (int i = 0; i <= 1024; i++) {
                Complex<Pow2.N4> r = Complex<Pow2.N4>.Ldexp((nu_abs + i) * v, 1);

                (a0, b0, a1, b1) = (a1, b1, r * a1 + a0, r * b1 + b0);
                s = a1 / b1;

                long exp = long.Max(a1.Exponent, b1.Exponent);
                if (exp <= int.MinValue) {
                    return Complex<Pow2.N4>.NaN;
                }

                (a0, a1, b0, b1) = (Complex<Pow2.N4>.Ldexp(a0, -exp), Complex<Pow2.N4>.Ldexp(a1, -exp), Complex<Pow2.N4>.Ldexp(b0, -exp), Complex<Pow2.N4>.Ldexp(b1, -exp));


                if (i > 0 && (i & 3) == 0) {
                    Complex<Pow2.N4> r0 = a0 * b1, r1 = a1 * b0;
                    if (!((r0 - r1).Magnitude > MultiPrecision<Pow2.N4>.Min(r0.Magnitude, r1.Magnitude) * 1e-36)) {
                        break;
                    }
                }
            }

            if (!Complex<Pow2.N4>.IsFinite(s)) {
                return 0d;
            }

            long exp_sum = 0;
            bool s_overone = s.Magnitude > 1d;
            (Complex<Pow2.N4> i0, Complex<Pow2.N4> i1) = s_overone ? (Complex<Pow2.N4>.One, 1d / s) : (s, 1d);

            for (int k = n - 1; k >= BesselUtil<Pow2.N4>.DirectMaxN - 1; k--) {
                (i1, i0) = (MultiPrecision<Pow2.N4>.Ldexp(k + alpha, 1) * v * i1 + i0, i1);

                long j0_exp = i0.Exponent, j1_exp = i1.Exponent;
                if (long.Sign(j0_exp) * long.Sign(j1_exp) > 0) {
                    long exp = j0_exp > 0 ? long.Max(j0_exp, j1_exp) : long.Min(j0_exp, j1_exp);
                    exp_sum += exp;
                    (i0, i1) = (Complex<Pow2.N4>.Ldexp(i0, -exp), Complex<Pow2.N4>.Ldexp(i1, -exp));
                }
            }

            Complex<Pow2.N4> y = Complex<Pow2.N4>.Ldexp(
                i0.Magnitude >= i1.Magnitude
                ? BesselN4.BesselI(alpha + (BesselUtil<Pow2.N4>.DirectMaxN - 1), z) / i0
                : BesselN4.BesselI(alpha + (BesselUtil<Pow2.N4>.DirectMaxN - 2), z) / i1,
                long.Clamp(-exp_sum, int.MinValue, int.MaxValue)
            ) * (s_overone ? 1d : s);

            if (MultiPrecision<Pow2.N4>.IsNegative(nu) && !MultiPrecision<Pow2.N4>.IsInteger(nu_abs)) {
                Complex<Pow2.N4> bk = 2d * MultiPrecision<Pow2.N4>.RcpPI * SinCosPICache<Pow2.N4>.SinPI(nu_abs) * BesselN4.BesselK(nu_abs, z);

                y += bk;
            }

            return y;
        }

        public static Complex<Pow2.N4> BesselK(MultiPrecision<Pow2.N4> nu, Complex<Pow2.N4> z) {
            nu = MultiPrecision<Pow2.N4>.Abs(nu);

            Debug.Assert(nu >= BesselUtil<Pow2.N4>.DirectMaxN);

            int n = (int)MultiPrecision<Pow2.N4>.Floor(nu);
            MultiPrecision<Pow2.N4> alpha = nu - n;

            Complex<Pow2.N4> k0 = BesselN4.BesselK(alpha + (BesselUtil<Pow2.N4>.DirectMaxN - 2), z);
            Complex<Pow2.N4> k1 = BesselN4.BesselK(alpha + (BesselUtil<Pow2.N4>.DirectMaxN - 1), z);

            if (Complex<Pow2.N4>.IsZero(k0) && Complex<Pow2.N4>.IsZero(k1)) {
                return 0d;
            }

            long exp_sum = 0;

            Complex<Pow2.N4> v = 1d / z;

            for (int k = BesselUtil<Pow2.N4>.DirectMaxN - 1; k < n; k++) {
                (k1, k0) = (MultiPrecision<Pow2.N4>.Ldexp(k + alpha, 1) * v * k1 + k0, k1);

                if (k1.Exponent > 0) {
                    long exp = k1.Exponent;
                    exp_sum += exp;
                    (k0, k1) = (Complex<Pow2.N4>.Ldexp(k0, -exp), Complex<Pow2.N4>.Ldexp(k1, -exp));
                }
            }

            k1 = Complex<Pow2.N4>.Ldexp(k1, long.Max(exp_sum, int.MinValue));

            return k1;
        }
    }
}
