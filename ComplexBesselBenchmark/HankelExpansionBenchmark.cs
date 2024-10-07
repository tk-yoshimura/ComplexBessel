using ComplexBessel;
using MultiPrecision;
using MultiPrecisionComplex;

namespace ComplexBesselBenchmark {
    [TestClass]
    public class HankelExpansionBenchmark {
        [TestMethod()]
        public void BesselJTest() {
            for (double nu = -16; nu <= 16; nu += 0.25) {
                using StreamWriter sw = new($"../../../../results/besselj_nu{nu}_hankel_convergence.csv");
                sw.WriteLine("r,i,relerr");

                for (double r = 0; r <= 64; r += 0.5) {
                    for (double i = 0; i <= 64; i += 0.5) {
                        if (r == 0 && i == 0) {
                            continue;
                        }

                        Complex<Pow2.N4> z4 = Limit<Pow2.N4>.BesselJ(nu, (r, i));

                        if (Complex<Pow2.N4>.IsNaN(z4)) {
                            sw.WriteLine($"{r},{i},1");
                            continue;
                        }

                        Complex<Pow2.N4> z5 = Limit<Plus1<Pow2.N4>>.BesselJ(nu, (r, i)).Convert<Pow2.N4>();

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

        [TestMethod()]
        public void BesselYTest() {
            for (double nu = -16; nu <= 16; nu += 0.25) {
                using StreamWriter sw = new($"../../../../results/bessely_nu{nu}_hankel_convergence.csv");
                sw.WriteLine("r,i,relerr");

                for (double r = 0; r <= 64; r += 0.5) {
                    for (double i = 0; i <= 64; i += 0.5) {
                        if (r == 0 && i == 0) {
                            continue;
                        }

                        Complex<Pow2.N4> z4 = Limit<Pow2.N4>.BesselY(nu, (r, i));

                        if (Complex<Pow2.N4>.IsNaN(z4)) {
                            sw.WriteLine($"{r},{i},1");
                            continue;
                        }

                        Complex<Pow2.N4> z5 = Limit<Plus1<Pow2.N4>>.BesselY(nu, (r, i)).Convert<Pow2.N4>();

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

        [TestMethod()]
        public void BesselITest() {
            for (double nu = -16; nu <= 16; nu += 0.25) {
                using StreamWriter sw = new($"../../../../results/besseli_nu{nu}_hankel_convergence.csv");
                sw.WriteLine("r,i,relerr");

                for (double r = 0; r <= 64; r += 0.5) {
                    for (double i = 0; i <= 64; i += 0.5) {
                        if (r == 0 && i == 0) {
                            continue;
                        }

                        Complex<Pow2.N4> z4 = Limit<Pow2.N4>.BesselI(nu, (r, i));

                        if (Complex<Pow2.N4>.IsNaN(z4)) {
                            sw.WriteLine($"{r},{i},1");
                            continue;
                        }

                        Complex<Pow2.N4> z5 = Limit<Plus1<Pow2.N4>>.BesselI(nu, (r, i)).Convert<Pow2.N4>();

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

        [TestMethod()]
        public void BesselKTest() {
            for (double nu = 0; nu <= 16; nu += 0.25) {
                using StreamWriter sw = new($"../../../../results/besselk_nu{nu}_hankel_convergence.csv");
                sw.WriteLine("r,i,relerr");

                for (double r = 0; r <= 64; r += 0.5) {
                    for (double i = 0; i <= 64; i += 0.5) {
                        if (r == 0 && i == 0) {
                            continue;
                        }

                        Complex<Pow2.N4> z4 = Limit<Pow2.N4>.BesselK(nu, (r, i));

                        if (Complex<Pow2.N4>.IsNaN(z4)) {
                            sw.WriteLine($"{r},{i},1");
                            continue;
                        }

                        Complex<Pow2.N4> z5 = Limit<Plus1<Pow2.N4>>.BesselK(nu, (r, i)).Convert<Pow2.N4>();

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