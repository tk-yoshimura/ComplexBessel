using DoubleDouble;

namespace DDoubleComplexBessel {
    static class SinCosPICache {
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
}