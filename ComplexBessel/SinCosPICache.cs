using MultiPrecision;

namespace ComplexBessel {
    static class SinCosPiCache<N> where N : struct, IConstant {
        private static readonly Dictionary<MultiPrecision<N>, MultiPrecision<N>> cospi_table = [];
        private static readonly Dictionary<MultiPrecision<N>, MultiPrecision<N>> sinpi_table = [];

        public static MultiPrecision<N> CosPi(MultiPrecision<N> theta) {
            if (!cospi_table.TryGetValue(theta, out MultiPrecision<N> cospi)) {
                cospi = MultiPrecision<N>.CosPi(theta);
                cospi_table[theta] = cospi;
            }

            return cospi;
        }

        public static MultiPrecision<N> SinPi(MultiPrecision<N> theta) {
            if (!sinpi_table.TryGetValue(theta, out MultiPrecision<N> sinpi)) {
                sinpi = MultiPrecision<N>.SinPi(theta);
                sinpi_table[theta] = sinpi;
            }

            return sinpi;
        }
    }
}