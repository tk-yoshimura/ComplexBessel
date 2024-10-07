using ComplexBessel;
using MultiPrecision;
using MultiPrecisionComplex;

namespace ComplexBesselBenchmark {
    [TestClass]
    public class HankelExpansionRange {
        [TestMethod()]
        public void BesselJTest() {
            using StreamWriter sw = new($"../../../../results/besselj_hankel_range.csv");

            sw.WriteLine("nu,norm");

            double max_norm = 0;

            for (double nu = -16; nu <= 16; nu += 0.25) {
                for (double norm = 32; norm <= 64; norm += 0.125) {
                    bool is_convergence = true;

                    for (double theta = 0; theta <= 0.5; theta += 1d / 32) {
                        Complex<Pow2.N4> z = (norm * MultiPrecision<Pow2.N4>.CosPI(theta), norm * MultiPrecision<Pow2.N4>.SinPI(theta));

                        Complex<Pow2.N4> z4 = Limit<Pow2.N4>.BesselJ(nu, z);

                        if (Complex<Pow2.N4>.IsNaN(z4)) {
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
            using StreamWriter sw = new($"../../../../results/bessely_hankel_range.csv");

            sw.WriteLine("nu,norm");

            double max_norm = 0;

            for (double nu = -16; nu <= 16; nu += 0.25) {
                for (double norm = 32; norm <= 64; norm += 0.125) {
                    bool is_convergence = true;

                    for (double theta = 0; theta <= 0.5; theta += 1d / 32) {
                        Complex<Pow2.N4> z = (norm * MultiPrecision<Pow2.N4>.CosPI(theta), norm * MultiPrecision<Pow2.N4>.SinPI(theta));

                        Complex<Pow2.N4> z4 = Limit<Pow2.N4>.BesselY(nu, z);

                        if (Complex<Pow2.N4>.IsNaN(z4)) {
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
            using StreamWriter sw = new($"../../../../results/besseli_hankel_range.csv");

            sw.WriteLine("nu,norm");

            double max_norm = 0;

            for (double nu = -16; nu <= 16; nu += 0.25) {
                for (double norm = 32; norm <= 64; norm += 0.125) {
                    bool is_convergence = true;

                    for (double theta = 0; theta <= 0.5; theta += 1d / 32) {
                        Complex<Pow2.N4> z = (norm * MultiPrecision<Pow2.N4>.CosPI(theta), norm * MultiPrecision<Pow2.N4>.SinPI(theta));

                        Complex<Pow2.N4> z4 = Limit<Pow2.N4>.BesselI(nu, z);

                        if (Complex<Pow2.N4>.IsNaN(z4)) {
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
            using StreamWriter sw = new($"../../../../results/besselk_hankel_range.csv");

            sw.WriteLine("nu,norm");

            double max_norm = 0;

            for (double nu = 0; nu <= 16; nu += 0.25) {
                for (double norm = 32; norm <= 64; norm += 0.125) {
                    bool is_convergence = true;

                    for (double theta = 0; theta <= 0.5; theta += 1d / 32) {
                        Complex<Pow2.N4> z = (norm * MultiPrecision<Pow2.N4>.CosPI(theta), norm * MultiPrecision<Pow2.N4>.SinPI(theta));

                        Complex<Pow2.N4> z4 = Limit<Pow2.N4>.BesselK(nu, z);

                        if (Complex<Pow2.N4>.IsNaN(z4)) {
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