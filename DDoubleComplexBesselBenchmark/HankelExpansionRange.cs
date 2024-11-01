using DDoubleComplexBessel;
using DoubleDouble;
using DoubleDoubleComplex;

namespace DDoubleComplexBesselBenchmark {
    [TestClass]
    public class HankelExpansionRange {
        private const double dnorm = 0.125;

        [TestMethod()]
        public void BesselJTest() {
            using StreamWriter sw = new($"../../../../results_ddouble/besselj_hankel_range.csv");

            sw.WriteLine("nu,norm");

            double max_norm = 0;

            for (double nu = -16; nu <= 16; nu += 0.25) {
                for (double norm = 50; norm >= 32; norm -= dnorm) {
                    bool is_convergence = true;

                    for (double theta = 0; theta <= 0.5; theta += 1d / 64) {
                        Complex z = (norm * ddouble.CosPi(theta), norm * ddouble.SinPi(theta));

                        Complex z4 = Limit.BesselJ(nu, z);

                        if (Complex.IsNaN(z4)) {
                            is_convergence = false;
                            break;
                        }
                    }

                    if (!is_convergence) {
                        max_norm = double.Max(norm + dnorm, max_norm);
                        sw.WriteLine($"{nu},{norm}");
                        break;
                    }
                }
            }

            sw.Close();

            Console.WriteLine(max_norm);
        }

        [TestMethod()]
        public void BesselYTest() {
            using StreamWriter sw = new($"../../../../results_ddouble/bessely_hankel_range.csv");

            sw.WriteLine("nu,norm");

            double max_norm = 0;

            for (double nu = -16; nu <= 16; nu += 0.25) {
                for (double norm = 50; norm >= 32; norm -= dnorm) {
                    bool is_convergence = true;

                    for (double theta = 0; theta <= 0.5; theta += 1d / 64) {
                        Complex z = (norm * ddouble.CosPi(theta), norm * ddouble.SinPi(theta));

                        Complex z4 = Limit.BesselY(nu, z);

                        if (Complex.IsNaN(z4)) {
                            is_convergence = false;
                            break;
                        }
                    }

                    if (!is_convergence) {
                        max_norm = double.Max(norm + dnorm, max_norm);
                        sw.WriteLine($"{nu},{norm}");
                        break;
                    }
                }
            }

            sw.Close();

            Console.WriteLine(max_norm);
        }

        [TestMethod()]
        public void BesselITest() {
            using StreamWriter sw = new($"../../../../results_ddouble/besseli_hankel_range.csv");

            sw.WriteLine("nu,norm");

            double max_norm = 0;

            for (double nu = -16; nu <= 16; nu += 0.25) {
                for (double norm = 50; norm >= 32; norm -= dnorm) {
                    bool is_convergence = true;

                    for (double theta = 0; theta <= 0.5; theta += 1d / 64) {
                        Complex z = (norm * ddouble.CosPi(theta), norm * ddouble.SinPi(theta));

                        Complex z4 = Limit.BesselI(nu, z);

                        if (Complex.IsNaN(z4)) {
                            is_convergence = false;
                            break;
                        }
                    }

                    if (!is_convergence) {
                        max_norm = double.Max(norm + dnorm, max_norm);
                        sw.WriteLine($"{nu},{norm}");
                        break;
                    }
                }
            }

            sw.Close();

            Console.WriteLine(max_norm);
        }

        [TestMethod()]
        public void BesselKTest() {
            using StreamWriter sw = new($"../../../../results_ddouble/besselk_hankel_range.csv");

            sw.WriteLine("nu,norm");

            double max_norm = 0;

            for (double nu = 0; nu <= 16; nu += 0.25) {
                for (double norm = 50; norm >= 32; norm -= dnorm) {
                    bool is_convergence = true;

                    for (double theta = 0; theta <= 0.5; theta += 1d / 64) {
                        Complex z = (norm * ddouble.CosPi(theta), norm * ddouble.SinPi(theta));

                        Complex z4 = Limit.BesselK(nu, z);

                        if (Complex.IsNaN(z4)) {
                            is_convergence = false;
                            break;
                        }
                    }

                    if (!is_convergence) {
                        max_norm = double.Max(norm + dnorm, max_norm);
                        sw.WriteLine($"{nu},{norm}");
                        break;
                    }
                }
            }

            sw.Close();

            Console.WriteLine(max_norm);
        }
    }
}