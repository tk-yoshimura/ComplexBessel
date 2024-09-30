using MultiPrecision;

namespace ComplexBessel {
    static class SinCosPICache<N> where N : struct, IConstant {
        private static readonly Dictionary<MultiPrecision<N>, MultiPrecision<N>> cospi_table = [];
        private static readonly Dictionary<MultiPrecision<N>, MultiPrecision<N>> sinpi_table = [];

        public static MultiPrecision<N> CosPI(MultiPrecision<N> theta) {
            if (!cospi_table.TryGetValue(theta, out MultiPrecision<N> cospi)) {
                cospi = MultiPrecision<N>.CosPI(theta);
                cospi_table[theta] = cospi;
            }

            return cospi;
        }

        public static MultiPrecision<N> SinPI(MultiPrecision<N> theta) {
            if (!sinpi_table.TryGetValue(theta, out MultiPrecision<N> sinpi)) {
                sinpi = MultiPrecision<N>.SinPI(theta);
                sinpi_table[theta] = sinpi;
            }

            return sinpi;
        }
    }
}