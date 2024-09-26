using ComplexBessel;
using MultiPrecision;
using MultiPrecisionComplex;

namespace ComplexBesselBenchmark {
    [TestClass]
    public class YoshidaPadeBenchmark {
        [TestMethod()]
        public void BesselKTest() {
            for (double nu = 0; nu <= 16; nu += 0.25) {
                YoshidaPade<Pow2.N4> pade_n4 = new(nu);
                YoshidaPade<Plus1<Pow2.N4>> pade_n5 = new(nu);

                using StreamWriter sw = new($"../../../../results/besselk_nu{nu}_yoshidapade_convergence.csv");
                sw.WriteLine("r,i,relerr");

                for (double r = 0; r <= 64; r += 0.5) {
                    for (double i = 0; i <= 64; i += 0.5) {
                        if (r == 0 && i == 0) {
                            continue;
                        }

                        Complex<Pow2.N4> z4 = pade_n4.BesselK((r, i));

                        if (Complex<Pow2.N4>.IsNaN(z4)) {
                            sw.WriteLine($"{r},{i},1");
                            continue;
                        }

                        Complex<Pow2.N4> z5 = pade_n5.BesselK((r, i)).Convert<Pow2.N4>();

                        if (Complex<Pow2.N4>.IsNaN(z5)) {
                            sw.WriteLine($"{r},{i},1");
                            continue;
                        }

                        MultiPrecision<Pow2.N4> err = (z4 - z5).Magnitude / z5.Magnitude;

                        sw.WriteLine($"{r},{i},{err:e10}");
                    }
                }
            }
        }
    }
}