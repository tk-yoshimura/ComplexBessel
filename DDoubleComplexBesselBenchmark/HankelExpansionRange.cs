using DDoubleComplexBessel;
using DoubleDouble;
using DoubleDoubleComplex;

namespace DDoubleComplexBesselBenchmark {
    [TestClass]
    public class HankelExpansionRange {
        [TestMethod()]
        public void BesselJTest() {
            using StreamWriter sw = new($"../../../../results_ddouble/besselj_hankel_range.csv");

            sw.WriteLine("nu,norm");

            double max_norm = 0;

            for (double nu = -16; nu <= 16; nu += 0.25) {
                HankelExpansion hankel_n4 = new(nu);

                for (double norm = 32; norm <= 64; norm += 0.125) {
                    bool is_convergence = true;

                    for (double theta = 0; theta <= 0.5; theta += 1d / 32) {
                        Complex z = (norm * ddouble.CosPI(theta), norm * ddouble.SinPI(theta));

                        Complex z4 = hankel_n4.BesselJ(z);

                        if (Complex.IsNaN(z4)) {
                            is_convergence = false;
                            break;
                        }
                    }

                    if (is_convergence) {
                        max_norm = double.Max(norm, max_norm);
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
                HankelExpansion hankel_n4 = new(nu);

                for (double norm = 32; norm <= 64; norm += 0.125) {
                    bool is_convergence = true;

                    for (double theta = 0; theta <= 0.5; theta += 1d / 32) {
                        Complex z = (norm * ddouble.CosPI(theta), norm * ddouble.SinPI(theta));

                        Complex z4 = hankel_n4.BesselY(z);

                        if (Complex.IsNaN(z4)) {
                            is_convergence = false;
                            break;
                        }
                    }

                    if (is_convergence) {
                        max_norm = double.Max(norm, max_norm);
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
                HankelExpansion hankel_n4 = new(nu);

                for (double norm = 32; norm <= 64; norm += 0.125) {
                    bool is_convergence = true;

                    for (double theta = 0; theta <= 0.5; theta += 1d / 32) {
                        Complex z = (norm * ddouble.CosPI(theta), norm * ddouble.SinPI(theta));

                        Complex z4 = hankel_n4.BesselI(z);

                        if (Complex.IsNaN(z4)) {
                            is_convergence = false;
                            break;
                        }
                    }

                    if (is_convergence) {
                        max_norm = double.Max(norm, max_norm);
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
                HankelExpansion hankel_n4 = new(nu);

                for (double norm = 32; norm <= 64; norm += 0.125) {
                    bool is_convergence = true;

                    for (double theta = 0; theta <= 0.5; theta += 1d / 32) {
                        Complex z = (norm * ddouble.CosPI(theta), norm * ddouble.SinPI(theta));

                        Complex z4 = hankel_n4.BesselK(z);

                        if (Complex.IsNaN(z4)) {
                            is_convergence = false;
                            break;
                        }
                    }

                    if (is_convergence) {
                        max_norm = double.Max(norm, max_norm);
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