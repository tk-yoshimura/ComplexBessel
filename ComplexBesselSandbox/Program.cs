using MultiPrecision;
using MultiPrecisionComplex;

namespace ComplexBesselSandbox {
    internal class Program {
        static void Main() {
            Complex<Pow2.N4> c = BesselY<Pow2.N4>(2, 2);

            Console.WriteLine("END");
            Console.Read();
        }

        public static Complex<N> BesselY<N>(int n, Complex<N> z) where N : struct, IConstant {
            Complex<N> c = 0;
            Complex<N> u = 1, v = 1, w = z * z / 4;

            MultiPrecision<N>[] fs = (new MultiPrecision<N>[n]).Select(
                (_, k) => MultiPrecision<N>.Gamma(n - k) / MultiPrecision<N>.Gamma(k + 1)
            ).ToArray();

            for (int k = 0; k < n; k++) {
                c += v * fs[k];
                v *= w;
            }
            c /= -v;

            Complex<N> h = 2 * (Complex<N>.Log(z / 2) + MultiPrecision<N>.EulerGamma);
            MultiPrecision<N> frac = 1 / MultiPrecision<N>.Gamma(n + 1);

            for (int k = 0, conv_times = 0; k < 256 && conv_times < 2; k++) {
                Complex<N> dc = u * frac * (h - MultiPrecision<N>.HarmonicNumber(k) - MultiPrecision<N>.HarmonicNumber(k + n));

                Complex<N> c_next = c + dc;

                if (c == c_next || !Complex<N>.IsFinite(c_next)) {
                    conv_times++;
                }
                else {
                    conv_times = 0;
                }

                c = c_next;
                u *= -w;
                frac /= (k + 1) * (n + k + 1);
            }

            Complex<N> y = c * MultiPrecision<N>.RcpPI * Complex<N>.Pow(z / 2, n);

            return y;
        }
    }
}
