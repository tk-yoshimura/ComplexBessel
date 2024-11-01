using MultiPrecision;

namespace DDoubleBesselYEpsMillerBackwardTests {
    public static class MPCoef<N> where N : struct, IConstant {
        public static (MultiPrecision<N> eta0, MultiPrecision<N> xi1) BesselYEta0Xi1(MultiPrecision<N> alpha, MultiPrecision<N> x) {
            MultiPrecision<N> s = MultiPrecision<N>.Pow(MultiPrecision<N>.Ldexp(1 / x, 1), alpha), sqs = s * s;

            MultiPrecision<N> rcot = 1d / MultiPrecision<N>.TanPi(alpha), rgamma = MultiPrecision<N>.Gamma(1d + alpha), rsqgamma = rgamma * rgamma;
            MultiPrecision<N> r = MultiPrecision<N>.Ldexp(MultiPrecision<N>.RcpPi * sqs, 1);
            MultiPrecision<N> p = sqs * rsqgamma * MultiPrecision<N>.RcpPi;

            MultiPrecision<N> eta0 = rcot - p / alpha;
            MultiPrecision<N> xi1 = rcot + p * (alpha * (alpha + 1d) + 1d) / (alpha * (alpha - 1d));

            return (eta0, xi1);
        }
    }
}
