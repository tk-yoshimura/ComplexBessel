using ComplexBessel;
using DDoubleComplexBessel;
using DoubleDouble;
using DoubleDoubleComplex;

namespace DDoubleComplexBesselBenchmark {
    [TestClass()]
    public class BesselTests {
        [TestMethod()]
        public void BesselJTest() {
            for (double nu = -16; nu <= 16; nu += 0.25) {
                using StreamWriter sw = new($"../../../../results_ddouble/besselj_nu{nu}_accuracy.csv");
                sw.WriteLine("r,i,z,relerr");

                for (double r = 0; r <= 42; r += 0.5) {
                    for (double i = 0; i <= 42; i += 0.5) {
                        if (r == 0 && i == 0) {
                            continue;
                        }

                        Complex expected = BesselN4.BesselJ(nu, (r, i)).ToString();
                        Complex actual = Bessel.BesselJ(nu, (r, i));

                        ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                        sw.WriteLine($"{r},{i},{actual},{err:e4}");
                    }
                }
            }
        }

        [TestMethod()]
        public void BesselYTest() {
            for (double nu = -16; nu <= 16; nu += 0.25) {
                using StreamWriter sw = new($"../../../../results_ddouble/bessely_nu{nu}_accuracy.csv");
                sw.WriteLine("r,i,z,relerr");

                for (double r = 0; r <= 42; r += 0.5) {
                    for (double i = 0; i <= 42; i += 0.5) {
                        if (r == 0 && i == 0) {
                            continue;
                        }

                        Complex expected = BesselN4.BesselY(nu, (r, i)).ToString();
                        Complex actual = Bessel.BesselY(nu, (r, i));

                        ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                        sw.WriteLine($"{r},{i},{actual},{err:e4}");
                    }
                }
            }

            for (int n = -16; n <= 16; n++) {
                foreach (double nu in new double[] { n - 0.125, n - 0.0625, n + 0.0625, n + 0.125 }) {
                    if (nu < -16 || nu > 16) {
                        continue;
                    }

                    using StreamWriter sw = new($"../../../../results_ddouble/bessely_nu{nu}_accuracy.csv");
                    sw.WriteLine("r,i,z,relerr");

                    for (double r = 0; r <= 42; r += 0.5) {
                        for (double i = 0; i <= 42; i += 0.5) {
                            if (r == 0 && i == 0) {
                                continue;
                            }

                            Complex expected = BesselN4.BesselY(nu, (r, i)).ToString();
                            Complex actual = Bessel.BesselY(nu, (r, i));

                            ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                            sw.WriteLine($"{r},{i},{actual},{err:e4}");
                        }
                    }
                }
            }
        }

        [TestMethod()]
        public void BesselITest() {
            for (double nu = -16; nu <= 16; nu += 0.25) {
                using StreamWriter sw = new($"../../../../results_ddouble/besseli_nu{nu}_accuracy.csv");
                sw.WriteLine("r,i,z,relerr");

                for (double r = 0; r <= 42; r += 0.5) {
                    for (double i = 0; i <= 42; i += 0.5) {
                        if (r == 0 && i == 0) {
                            continue;
                        }

                        Complex expected = BesselN4.BesselI(nu, (r, i)).ToString();
                        Complex actual = Bessel.BesselI(nu, (r, i));

                        ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                        sw.WriteLine($"{r},{i},{actual},{err:e4}");
                    }
                }
            }
        }

        [TestMethod()]
        public void BesselKTest() {
            for (double nu = 0; nu <= 16; nu += 0.25) {
                using StreamWriter sw = new($"../../../../results_ddouble/besselk_nu{nu}_accuracy.csv");
                sw.WriteLine("r,i,z,relerr");

                for (double r = 0; r <= 42; r += 0.5) {
                    for (double i = 0; i <= 42; i += 0.5) {
                        if (r == 0 && i == 0) {
                            continue;
                        }

                        Complex expected = BesselN4.BesselK(nu, (r, i)).ToString();
                        Complex actual = Bessel.BesselK(nu, (r, i));

                        ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                        sw.WriteLine($"{r},{i},{actual},{err:e4}");
                    }
                }
            }

            for (int n = 0; n <= 16; n++) {
                foreach (double nu in new double[] { n - 0.125, n - 0.0625, n + 0.0625, n + 0.125 }) {
                    if (nu < 0 || nu > 16) {
                        continue;
                    }

                    using StreamWriter sw = new($"../../../../results_ddouble/besselk_nu{nu}_accuracy.csv");
                    sw.WriteLine("r,i,z,relerr");

                    for (double r = 0; r <= 42; r += 0.5) {
                        for (double i = 0; i <= 42; i += 0.5) {
                            if (r == 0 && i == 0) {
                                continue;
                            }

                            Complex expected = BesselN4.BesselK(nu, (r, i)).ToString();
                            Complex actual = Bessel.BesselK(nu, (r, i));

                            ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                            sw.WriteLine($"{r},{i},{actual},{err:e4}");
                        }
                    }
                }
            }
        }
    }
}