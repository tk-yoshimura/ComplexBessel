using DoubleDouble;

namespace DDoubleComplexBessel {
    static class SinCosPiCache {
        private static readonly Dictionary<ddouble, ddouble> cospi_table = [];
        private static readonly Dictionary<ddouble, ddouble> sinpi_table = [];

        public static ddouble CosPi(ddouble theta) {
            if (!cospi_table.TryGetValue(theta, out ddouble cospi)) {
                cospi = ddouble.CosPi(theta);
                cospi_table[theta] = cospi;
            }

            return cospi;
        }

        public static ddouble SinPi(ddouble theta) {
            if (!sinpi_table.TryGetValue(theta, out ddouble sinpi)) {
                sinpi = ddouble.SinPi(theta);
                sinpi_table[theta] = sinpi;
            }

            return sinpi;
        }
    }
}