using ComplexBessel;
using MultiPrecision;
using MultiPrecisionComplex;

namespace ComplexBesselBenchmark {
    [TestClass]
    public class PowerSeriesRange {
        [TestMethod()]
        public void BesselJTest() {
            using StreamWriter sw = new($"../../../../results/besselj_powerseries_range.csv");
            sw.WriteLine("r,i");

            double r = 0;

            for (double i = 0; i <= 64; i += 0.125) {
                for (r = double.Max(0, r - 4); r <= 64; r += 0.125) {
                    if (r == 0 && i == 0) {
                        continue;
                    }

                    Complex<Pow2.N4> z4 = PowerSeries<Pow2.N4>.BesselJ(0, (r, i));

                    if (Complex<Pow2.N4>.IsNaN(z4)) {
                        sw.WriteLine($"{r},{i},1");
                        continue;
                    }

                    Complex<Pow2.N4> z5 = PowerSeries<Plus1<Pow2.N4>>.BesselJ(0, (r, i)).Convert<Pow2.N4>();

                    if (Complex<Pow2.N4>.IsNaN(z5)) {
                        sw.WriteLine($"{r},{i},1");
                        continue;
                    }

                    MultiPrecision<Pow2.N4> err = (z4 - z5).Magnitude / z5.Magnitude;

                    if (err > 3.16e-37) { 
                        sw.WriteLine($"{r - 0.125},{i}");
                        break;
                    }
                }
            }
        }
    }
}