using ComplexBessel;
using DDoubleComplexBessel;
using DoubleDouble;
using DoubleDoubleComplex;

namespace DDoubleComplexBesselBenchmark {
    [TestClass]
    public class PowerSeriesRange {
        [TestMethod()]
        public void BesselJTest() {
            for (double nu = 0; nu <= 16; nu += 0.25) {
                using StreamWriter sw = new($"../../../../results_ddouble/besselj_nu{nu}_powerseries_range.csv");
                sw.WriteLine("r,i");

                double r = 0;

                for (double i = 0; i <= 64; i += 0.125) {
                    for (r = double.Max(0, r - 4); r <= 64; r += 0.125) {
                        if (r == 0 && i == 0) {
                            continue;
                        }

                        Complex z4 = PowerSeries.BesselJ(nu, (r, i));

                        if (Complex.IsNaN(z4)) {
                            sw.WriteLine($"{r},{i},1");
                            continue;
                        }

                        Complex z5 = BesselN4.BesselJ(nu, (r, i)).ToString();

                        if (Complex.IsNaN(z5)) {
                            sw.WriteLine($"{r},{i},1");
                            continue;
                        }

                        ddouble err = (z4 - z5).Magnitude / z5.Magnitude;

                        if (err > 3.16e-30) {
                            sw.WriteLine($"{r - 0.125},{i}");
                            break;
                        }
                    }

                    if (r * r + i * i >= 40 * 40) {
                        break;
                    }
                }
            }
        }
    }
}