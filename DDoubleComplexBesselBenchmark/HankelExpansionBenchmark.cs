using ComplexBessel;
using DDoubleComplexBessel;
using DoubleDouble;
using DoubleDoubleComplex;

namespace DDoubleComplexBesselBenchmark {
    [TestClass]
    public class HankelExpansionBenchmark {
        [TestMethod()]
        public void BesselJTest() {
            for (double nu = -16; nu <= 16; nu += 0.25) {
                using StreamWriter sw = new($"../../../../results_ddouble/besselj_nu{nu}_hankel_convergence.csv");
                sw.WriteLine("r,i,relerr");

                for (double r = 0; r <= 42; r += 0.5) {
                    for (double i = 0; i <= 42; i += 0.5) {
                        if (r == 0 && i == 0) {
                            continue;
                        }

                        Complex zdd = Limit.BesselJ(nu, (r, i));

                        if (Complex.IsNaN(zdd)) {
                            sw.WriteLine($"{r},{i},1");
                            continue;
                        }

                        Complex zmp = BesselN4.BesselJ(nu, (r, i)).ToString();

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

        [TestMethod()]
        public void BesselYTest() {
            for (double nu = -16; nu <= 16; nu += 0.25) {
                using StreamWriter sw = new($"../../../../results_ddouble/bessely_nu{nu}_hankel_convergence.csv");
                sw.WriteLine("r,i,relerr");

                for (double r = 0; r <= 42; r += 0.5) {
                    for (double i = 0; i <= 42; i += 0.5) {
                        if (r == 0 && i == 0) {
                            continue;
                        }

                        Complex zdd = Limit.BesselY(nu, (r, i));

                        if (Complex.IsNaN(zdd)) {
                            sw.WriteLine($"{r},{i},1");
                            continue;
                        }

                        Complex zmp = BesselN4.BesselY(nu, (r, i)).ToString();

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

        [TestMethod()]
        public void BesselITest() {
            for (double nu = -16; nu <= 16; nu += 0.25) {
                using StreamWriter sw = new($"../../../../results_ddouble/besseli_nu{nu}_hankel_convergence.csv");
                sw.WriteLine("r,i,relerr");

                for (double r = 0; r <= 42; r += 0.5) {
                    for (double i = 0; i <= 42; i += 0.5) {
                        if (r == 0 && i == 0) {
                            continue;
                        }

                        Complex zdd = Limit.BesselI(nu, (r, i));

                        if (Complex.IsNaN(zdd)) {
                            sw.WriteLine($"{r},{i},1");
                            continue;
                        }

                        Complex zmp = BesselN4.BesselI(nu, (r, i)).ToString();

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

        [TestMethod()]
        public void BesselKTest() {
            for (double nu = 0; nu <= 16; nu += 0.25) {
                using StreamWriter sw = new($"../../../../results_ddouble/besselk_nu{nu}_hankel_convergence.csv");
                sw.WriteLine("r,i,relerr");

                for (double r = 0; r <= 42; r += 0.5) {
                    for (double i = 0; i <= 42; i += 0.5) {
                        if (r == 0 && i == 0) {
                            continue;
                        }

                        Complex zdd = Limit.BesselK(nu, (r, i));

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