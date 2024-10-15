using MultiPrecision;

namespace ComplexBessel {
    static class BesselUtil<N> where N : struct, IConstant {
        public static readonly double Eps = double.ScaleB(1, -MultiPrecision<N>.Bits);
        public static readonly double ExtremelyNearZero = double.ScaleB(1, -28);
        public static readonly double InterpolationThreshold = double.ScaleB(1, -25);
        public static readonly double MillerBwdBesselYEps = double.ScaleB(1, -30);

        public const int RecurrenceMaxN = 256;
        public const int DirectMaxN = 16;

        public static void CheckNu(MultiPrecision<N> nu) {
            if (!(MultiPrecision<N>.Abs(nu) <= RecurrenceMaxN)) {
                throw new ArgumentOutOfRangeException(
                    nameof(nu),
                    $"In the calculation of the Bessel function, nu with an absolute value greater than {RecurrenceMaxN} is not supported."
                );
            }
        }

        public static void CheckN(int n) {
            if (n < -RecurrenceMaxN || n > RecurrenceMaxN) {
                throw new ArgumentOutOfRangeException(
                    nameof(n),
                    $"In the calculation of the Bessel function, n with an absolute value greater than {RecurrenceMaxN} is not supported."
                );
            }
        }

        public static bool UseRecurrence(MultiPrecision<N> nu) {
            return MultiPrecision<N>.Abs(nu) > DirectMaxN;
        }

        public static bool NearlyInteger(MultiPrecision<N> nu, out int n) {
            n = (int)MultiPrecision<N>.Round(nu);

            return MultiPrecision<N>.Abs(nu - n) < Eps;
        }
    }
}
