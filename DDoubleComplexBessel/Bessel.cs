using DoubleDouble;
using DoubleDoubleComplex;

namespace DDoubleComplexBessel {
    public static class Bessel {
        private const double hankel_threshold = 38.875, miller_backward_threshold = 6;
        private const double besselk_pade_threshold = 2, besselk_nz_threshold = 2, besseljy_powerseries_bias = 2;

        public static Complex BesselJ(ddouble nu, Complex z) {
            BesselUtil.CheckNu(nu);

            if (ddouble.IsNegative(z.I)) {
                return BesselJ(nu, z.Conj).Conj;
            }

            if (ddouble.IsNegative(z.R)) {
                return (SinCosPiCache.CosPi(nu), SinCosPiCache.SinPi(nu)) * BesselJ(nu, (-z.R, z.I)).Conj;
            }

            if (z.Magnitude >= hankel_threshold) {
                return Limit.BesselJ(nu, z);
            }
            else if (z.R <= PowerSeriesThreshold(nu, z.I)) {
                return PowerSeries.BesselJ(nu, z);
            }
            else if (z.I <= miller_backward_threshold) {
                return MillerBackward.BesselJ(nu, z);
            }
            else {
                return (SinCosPiCache.CosPi(nu / 2), SinCosPiCache.SinPi(nu / 2)) * BesselI(nu, (z.I, z.R)).Conj;
            }
        }

        public static Complex BesselY(ddouble nu, Complex z) {
            BesselUtil.CheckNu(nu);

            if (ddouble.IsNegative(z.I)) {
                return BesselY(nu, z.Conj).Conj;
            }

            if (ddouble.IsNegative(z.R)) {
                return (SinCosPiCache.CosPi(nu), -SinCosPiCache.SinPi(nu)) * BesselY(nu, (-z.R, z.I)).Conj
                     + (0, 2 * SinCosPiCache.CosPi(nu)) * BesselJ(nu, (-z.R, z.I)).Conj;
            }

            if (z.Magnitude >= hankel_threshold) {
                return Limit.BesselY(nu, z);
            }
            else if (z.R <= PowerSeriesThreshold(nu, z.I) - besseljy_powerseries_bias) {
                if (ddouble.Round(nu) == nu || ddouble.ILogB(ddouble.Round(nu) - nu) >= -3) {
                    return PowerSeries.BesselY(nu, z);
                }
                else if (z.I <= miller_backward_threshold / 2) {
                    return MillerBackward.BesselY(nu, z);
                }
            }
            else if (z.I <= miller_backward_threshold) {
                return MillerBackward.BesselY(nu, z);
            }
            {
                Complex c = (SinCosPiCache.CosPi(nu / 2), SinCosPiCache.SinPi(nu / 2));
                Complex bi = BesselI(nu, (z.I, z.R));
                Complex bk = BesselK(nu, (z.I, z.R));

                Complex y = Complex.ImaginaryOne * c * bi.Conj - 2 * ddouble.RcpPi * (c * bk).Conj;

                return y;
            }
        }

        public static Complex BesselI(ddouble nu, Complex z) {
            BesselUtil.CheckNu(nu);

            if (ddouble.IsNegative(z.I)) {
                return BesselI(nu, z.Conj).Conj;
            }

            if (ddouble.IsNegative(z.R)) {
                return (SinCosPiCache.CosPi(nu), SinCosPiCache.SinPi(nu)) * BesselI(nu, (-z.R, z.I)).Conj;
            }

            if (z.Magnitude >= hankel_threshold) {
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
            BesselUtil.CheckNu(nu);

            nu = ddouble.Abs(nu);

            if (ddouble.IsNegative(z.I)) {
                return BesselK(nu, z.Conj).Conj;
            }

            if (ddouble.IsNegative(z.R)) {
                return (SinCosPiCache.CosPi(nu), -SinCosPiCache.SinPi(nu)) * BesselK(nu, (-z.R, z.I)).Conj
                     - (0, ddouble.Pi) * BesselI(nu, (-z.R, z.I)).Conj;
            }

            if (z.Magnitude >= hankel_threshold) {
                return Limit.BesselK(nu, z);
            }
            else if (z.Magnitude <= besselk_nz_threshold) {
                if (ddouble.Round(nu) == nu || ddouble.ILogB(ddouble.Round(nu) - nu) >= -3) {
                    return PowerSeries.BesselK(nu, z);
                }
                else {
                    return AmosPowerSeries.BesselK(nu, z);
                }
            }
            else if (z.R >= besselk_pade_threshold) {
                return YoshidaPade.BesselK(nu, z);
            }
            else {
                Complex c = (SinCosPiCache.CosPi(nu / 2), -SinCosPiCache.SinPi(nu / 2));
                Complex bi = BesselI(nu, z);
                Complex by = BesselY(nu, (z.I, z.R));

                Complex y = c * (-Complex.ImaginaryOne * c * bi - by.Conj) * ddouble.Pi / 2;

                return y;
            }
        }

        static ddouble PowerSeriesThreshold(ddouble nu, ddouble x) {
            ddouble nu_abs = ddouble.Abs(nu);

            return 7.5 + nu_abs * (3.57e-1 + nu_abs * 5.23e-3) + x * (4.67e-1 - nu_abs * 1.51e-2);
        }
    }
}
