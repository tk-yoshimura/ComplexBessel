using ComplexBessel;
using MultiPrecision;
using MultiPrecisionComplex;

namespace ComplexBesselBenchmark {
    [TestClass]
    public class AmosPowerSeriesBenchmark {
        [TestMethod()]
        public void BesselKTest() {
            for (double nu = 0; nu <= 0.5; nu = nu < double.ScaleB(1, -16) ? double.ScaleB(1, -16) : nu * 2) {
                using StreamWriter sw = new($"../../../../results/besselk_nu{nu}_amos_besselk_convergence.csv");
                sw.WriteLine("r,i,relerr");

                for (double r = 0; r <= 8; r += 0.125) {
                    for (double i = 0; i <= 8; i += 0.125) {
                        if (r == 0 && i == 0) {
                            continue;
                        }

                        Complex<Pow2.N4> z4 = AmosPowerSeries<Pow2.N4>.BesselK(nu, (r, i));

                        if (Complex<Pow2.N4>.IsNaN(z4)) {
                            sw.WriteLine($"{r},{i},1");
                            continue;
                        }

                        Complex<Pow2.N4> z5 = AmosPowerSeries<Plus1<Pow2.N4>>.BesselK(nu, (r, i)).Convert<Pow2.N4>();

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