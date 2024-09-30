using ComplexBessel;
using DDoubleComplexBessel;
using DoubleDouble;
using DoubleDoubleComplex;

namespace DDoubleComplexBesselBenchmark {
    [TestClass]
    public class YoshidaPadeBenchmark {
        [TestMethod()]
        public void BesselKTest() {
            for (double nu = 0; nu <= 16; nu += 0.25) {
                YoshidaPade pade_dd = new(nu);

                using StreamWriter sw = new($"../../../../results_ddouble/besselk_nu{nu}_yoshidapade_convergence.csv");
                sw.WriteLine("r,i,relerr");

                for (double r = 0; r <= 42; r += 0.5) {
                    for (double i = 0; i <= 42; i += 0.5) {
                        if (r == 0 && i == 0) {
                            continue;
                        }

                        Complex zdd = pade_dd.BesselK((r, i));

                        if (Complex.IsNaN(zdd)) {
                            sw.WriteLine($"{r},{i},1");
                            continue;
                        }

                        Complex zmp = BesselN4.BesselK(nu, (r, i)).ToString();

                        if (Complex.IsNaN(zmp)) {
                            sw.WriteLine($"{r},{i},1");
                            continue;
                        }

                        ddouble err = (zdd - zmp).Magnitude / zmp.Magnitude;

                        sw.WriteLine($"{r},{i},{err:e10}");
                    }
                }
            }
        }
    }
}