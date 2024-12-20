﻿using MultiPrecision;
using CMP4 = MultiPrecisionComplex.Complex<MultiPrecision.Pow2.N4>;
using MP4 = MultiPrecision.MultiPrecision<MultiPrecision.Pow2.N4>;

namespace ComplexBessel {
    public static class BesselN4 {
        public static CMP4 BesselJ(MP4 nu, CMP4 z) {
            BesselUtil<Pow2.N4>.CheckNu(nu);

            if (z.I.Sign == Sign.Minus) {
                return BesselJ(nu, z.Conj).Conj;
            }

            if (z.R.Sign == Sign.Minus) {
                return (SinCosPiCache<Pow2.N4>.CosPi(nu), SinCosPiCache<Pow2.N4>.SinPi(nu)) * BesselJ(nu, (-z.R, z.I)).Conj;
            }

            if (BesselUtil<Pow2.N4>.UseRecurrence(nu)) {
                return Recurrence.BesselJ(nu, z);
            }

            if (z.Magnitude >= 46.5) {
                return Limit<Pow2.N4>.BesselJ(nu, z);
            }
            else if (z.R <= PowerSeriesThreshold(nu, z.I)) {
                return PowerSeries<Pow2.N4>.BesselJ(nu, z);
            }
            else if (z.I <= 8) {
                return MillerBackward<Pow2.N4>.BesselJ(nu, z);
            }
            else {
                return (SinCosPiCache<Pow2.N4>.CosPi(nu / 2), SinCosPiCache<Pow2.N4>.SinPi(nu / 2)) * BesselI(nu, (z.I, z.R)).Conj;
            }
        }

        public static CMP4 BesselY(MP4 nu, CMP4 z) {
            BesselUtil<Pow2.N4>.CheckNu(nu);

            if (z.I.Sign == Sign.Minus) {
                return BesselY(nu, z.Conj).Conj;
            }

            if (z.R.Sign == Sign.Minus) {
                return (SinCosPiCache<Pow2.N4>.CosPi(nu), -SinCosPiCache<Pow2.N4>.SinPi(nu)) * BesselY(nu, (-z.R, z.I)).Conj + (0, 2 * SinCosPiCache<Pow2.N4>.CosPi(nu)) * BesselJ(nu, (-z.R, z.I)).Conj;
            }

            if (BesselUtil<Pow2.N4>.UseRecurrence(nu)) {
                return Recurrence.BesselY(nu, z);
            }

            if (z.Magnitude >= 46.5) {
                return Limit<Pow2.N4>.BesselY(nu, z);
            }
            else if (z.R <= PowerSeriesThreshold(nu, z.I)) {
                if (MP4.Round(nu) == nu || (MP4.Round(nu) - nu).Exponent >= -4) {
                    return PowerSeries<Pow2.N4>.BesselY(nu, z);
                }
                else if (z.I <= 4) {
                    return MillerBackward<Pow2.N4>.BesselY(nu, z);
                }
            }
            else if (z.I <= 8) {
                return MillerBackward<Pow2.N4>.BesselY(nu, z);
            }
            {
                CMP4 c = (SinCosPiCache<Pow2.N4>.CosPi(nu / 2), SinCosPiCache<Pow2.N4>.SinPi(nu / 2));
                CMP4 bi = BesselI(nu, (z.I, z.R));
                CMP4 bk = BesselK(nu, (z.I, z.R));

                CMP4 y = (0, 1) * c * bi.Conj - 2 * MP4.RcpPi * (c * bk).Conj;

                return y;
            }
        }

        public static CMP4 BesselI(MP4 nu, CMP4 z) {
            BesselUtil<Pow2.N4>.CheckNu(nu);

            if (z.I.Sign == Sign.Minus) {
                return BesselI(nu, z.Conj).Conj;
            }

            if (z.R.Sign == Sign.Minus) {
                return (SinCosPiCache<Pow2.N4>.CosPi(nu), SinCosPiCache<Pow2.N4>.SinPi(nu)) * BesselI(nu, (-z.R, z.I)).Conj;
            }

            if (BesselUtil<Pow2.N4>.UseRecurrence(nu)) {
                return Recurrence.BesselI(nu, z);
            }

            if (z.Magnitude >= 46.25) {
                return Limit<Pow2.N4>.BesselI(nu, z);
            }
            else if (z.I <= PowerSeriesThreshold(nu, z.R)) {
                return PowerSeries<Pow2.N4>.BesselI(nu, z);
            }
            else {
                return MillerBackward<Pow2.N4>.BesselI(nu, z);
            }
        }

        public static CMP4 BesselK(MP4 nu, CMP4 z) {
            BesselUtil<Pow2.N4>.CheckNu(nu);

            nu = MP4.Abs(nu);

            if (z.I.Sign == Sign.Minus) {
                return BesselK(nu, z.Conj).Conj;
            }

            if (z.R.Sign == Sign.Minus) {
                return (SinCosPiCache<Pow2.N4>.CosPi(nu), -SinCosPiCache<Pow2.N4>.SinPi(nu)) * BesselK(nu, (-z.R, z.I)).Conj - (0, MP4.Pi) * BesselI(nu, (-z.R, z.I)).Conj;
            }

            if (BesselUtil<Pow2.N4>.UseRecurrence(nu)) {
                return Recurrence.BesselK(nu, z);
            }

            if (z.Magnitude >= 44.5) {
                return Limit<Pow2.N4>.BesselK(nu, z);
            }
            else if (z.Magnitude <= 4) {
                if (MP4.Round(nu) == nu || (MP4.Round(nu) - nu).Exponent >= -4) {
                    return PowerSeries<Pow2.N4>.BesselK(nu, z);
                }
                else {
                    return AmosPowerSeries<Pow2.N4>.BesselK(nu, z);
                }
            }
            else if (z.R >= MP4.Min(2, z.I / 2)) {
                return YoshidaPade<Pow2.N4>.BesselK(nu, z);
            }
            else {
                CMP4 c = (SinCosPiCache<Pow2.N4>.CosPi(nu / 2), -SinCosPiCache<Pow2.N4>.SinPi(nu / 2));
                CMP4 bi = BesselI(nu, z);
                CMP4 by = BesselY(nu, (z.I, z.R));

                CMP4 y = c * ((0, -1) * c * bi - by.Conj) * MP4.Pi / 2;

                return y;
            }
        }

        public static CMP4 HankelH1(MP4 nu, CMP4 z) {
            if (MP4.IsNegative(z.R) && MP4.IsNegative(z.I)) {
                return BesselJ(nu, z) + (0, 1) * BesselY(nu, z);
            }
            else {
                CMP4 c = (-SinCosPiCache<Pow2.N4>.SinPi(nu / 2), -SinCosPiCache<Pow2.N4>.CosPi(nu / 2));
                CMP4 f = BesselK(nu, (z.I, z.R)).Conj;

                CMP4 y = 2 * MP4.RcpPi * c * f;

                return y;
            }
        }

        public static CMP4 HankelH2(MP4 nu, CMP4 z) {
            if (MP4.IsNegative(z.R) && MP4.IsPositive(z.I)) {
                return BesselJ(nu, z) - (0, 1) * BesselY(nu, z);
            }
            else {
                CMP4 c = (-SinCosPiCache<Pow2.N4>.SinPi(nu / 2), SinCosPiCache<Pow2.N4>.CosPi(nu / 2));
                CMP4 f = BesselK(nu, (-z.I, z.R));

                CMP4 y = 2 * MP4.RcpPi * c * f;

                return y;
            }
        }

        static MP4 PowerSeriesThreshold(MP4 nu, MP4 x) {
            MP4 nu_abs = MP4.Abs(nu);
            return 6 + nu_abs * 4.59e-1 + x * (4.20e-1 - nu_abs * 1.44e-2);
        }
    }
}
