using DoubleDouble;
using System.Collections.Concurrent;

namespace DDoubleComplexBessel {
    static class SinCosPiCache {
        private static readonly ConcurrentDictionary<ddouble, ddouble> cospi_table = [], sinpi_table = [];

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