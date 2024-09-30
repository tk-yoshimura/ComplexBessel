using DoubleDouble;

namespace DDoubleComplexBessel {
    static class BesselUtil {
        public static readonly double Eps = double.ScaleB(1, -1000);
        public static readonly double ExtremelyNearZero = double.ScaleB(1, -28);
        public static readonly double InterpolationThreshold = double.ScaleB(1, -25);
        public static readonly double MillerBwdBesselYEps = double.ScaleB(1, -30);

        public const int MaxN = 16;

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
    }
}
