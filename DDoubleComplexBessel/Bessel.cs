using DoubleDouble;
using DoubleDoubleComplex;

namespace DDoubleComplexBessel {
    public static class Bessel {
        private const double hankel_threshold = 38.75, miller_backward_threshold = 6;
        private const double besselk_pade_threshold = 1, besselk_nz_threshold = 4;
        private static readonly Dictionary<ddouble, HankelExpansion> hankel_table = [];
        private static readonly Dictionary<ddouble, YoshidaPade> pade_table = [];

        public static Complex BesselJ(ddouble nu, Complex z) {
            BesselUtil.CheckNu(nu);

            if (ddouble.IsNegative(z.I)) {
                return BesselJ(nu, z.Conj).Conj;
            }

            if (ddouble.IsNegative(z.R)) {
                return (SinCosPICache.CosPI(nu), SinCosPICache.SinPI(nu)) * BesselJ(nu, -z);
            }

            if (z.Magnitude >= hankel_threshold) {
                if (!hankel_table.TryGetValue(nu, out HankelExpansion hankel)) {
                    hankel = new(nu);
                    hankel_table[nu] = hankel;
                }

                return hankel.BesselJ(z);
            }
            else if (z.R <= PowerSeriesThreshold(nu, z.I)) {
                return PowerSeries.BesselJ(nu, z);
            }
            else if (z.I <= miller_backward_threshold) {
                return MillerBackward.BesselJ(nu, z);
            }
            else {
                return (SinCosPICache.CosPI(nu / 2), SinCosPICache.SinPI(nu / 2)) * BesselI(nu, (z.I, z.R)).Conj;
            }
        }

        public static Complex BesselY(ddouble nu, Complex z) {
            BesselUtil.CheckNu(nu);

            if (ddouble.IsNegative(z.I)) {
                return BesselY(nu, z.Conj).Conj;
            }

            if (ddouble.IsNegative(z.R)) {
                return (SinCosPICache.CosPI(nu), -SinCosPICache.SinPI(nu)) * BesselY(nu, -z) 
                     + (0, 2 * SinCosPICache.CosPI(nu)) * BesselJ(nu, -z);
            }

            if (z.Magnitude >= hankel_threshold) {
                if (!hankel_table.TryGetValue(nu, out HankelExpansion hankel)) {
                    hankel = new(nu);
                    hankel_table[nu] = hankel;
                }

                return hankel.BesselY(z);
            }
            else if (z.R <= PowerSeriesThreshold(nu, z.I) - 2) {
                return PowerSeries.BesselY(nu, z);
            }
            else if (z.I <= miller_backward_threshold) {
                return MillerBackward.BesselY(nu, z);
            }
            else {
                Complex c = (SinCosPICache.CosPI(nu / 2), SinCosPICache.SinPI(nu / 2));
                Complex bi = BesselI(nu, (z.I, z.R));
                Complex bk = BesselK(nu, (z.I, z.R));

                Complex y = (0, 1) * c * bi.Conj - 2 * ddouble.RcpPI * (c * bk).Conj;

                return y;
            }
        }

        public static Complex BesselI(ddouble nu, Complex z) {
            BesselUtil.CheckNu(nu);

            if (ddouble.IsNegative(z.I)) {
                return BesselI(nu, z.Conj).Conj;
            }

            if (ddouble.IsNegative(z.R)) {
                return (SinCosPICache.CosPI(nu), SinCosPICache.SinPI(nu)) * BesselI(nu, -z);
            }

            if (z.Magnitude >= hankel_threshold) {
                if (!hankel_table.TryGetValue(nu, out HankelExpansion hankel)) {
                    hankel = new(nu);
                    hankel_table[nu] = hankel;
                }

                return hankel.BesselI(z);
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
                return (SinCosPICache.CosPI(nu), -SinCosPICache.SinPI(nu)) * BesselK(nu, -z)
                     - (0, ddouble.PI) * BesselI(nu, -z);
            }

            if (z.Magnitude >= hankel_threshold) {
                if (!hankel_table.TryGetValue(nu, out HankelExpansion hankel)) {
                    hankel = new(nu);
                    hankel_table[nu] = hankel;
                }

                return hankel.BesselK(z);
            }
            else if (z.Magnitude <= besselk_nz_threshold) {
                return PowerSeries.BesselK(nu, z);
            }
            else if (z.R >= besselk_pade_threshold) {
                if (!pade_table.TryGetValue(nu, out YoshidaPade pade)) {
                    pade = new(nu);
                    pade_table[nu] = pade;
                }

                return pade.BesselK(z);
            }
            else {
                Complex c = (SinCosPICache.CosPI(nu / 2), -SinCosPICache.SinPI(nu / 2));
                Complex bi = BesselI(nu, z);
                Complex by = BesselY(nu, (z.I, z.R));

                Complex y = c * ((0, -1) * c * bi - by.Conj) * ddouble.PI / 2;

                return y;
            }
        }

        static ddouble PowerSeriesThreshold(ddouble nu, ddouble x) {
            ddouble nu_abs = ddouble.Abs(nu);

            return 7.5 + nu_abs * (3.57e-1 + nu_abs * 5.23e-3) + x * (4.67e-1 - nu_abs * 1.51e-2);
        }
    }
}
