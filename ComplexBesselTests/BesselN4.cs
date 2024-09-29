using ComplexBessel;
using MultiPrecision;
using MultiPrecisionComplex;

namespace ComplexBesselTests {
    [TestClass()]
    public class BesselN4Tests {
        [TestMethod()]
        public void BesselJTest() {
            for (double nu = -16; nu <= 16; nu += 0.25) {
                Console.WriteLine(nu);

                foreach (double r in new double[] { 0, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64 }) {
                    foreach (double i in new double[] { 0, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64 }) {
                        if (r == 0 && i == 0) {
                            continue;
                        }

                        Complex<Pow2.N4> z = (r, i);

                        Console.WriteLine($"{z}: {BesselN4.BesselJ(nu, z)}");
                    }
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselYTest() {
            for (double nu = -16; nu <= 16; nu += 0.25) {
                Console.WriteLine(nu);

                foreach (double r in new double[] { 0, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64 }) {
                    foreach (double i in new double[] { 0, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64 }) {
                        if (r == 0 && i == 0) {
                            continue;
                        }

                        Complex<Pow2.N4> z = (r, i);

                        Console.WriteLine($"{z}: {BesselN4.BesselY(nu, z)}");
                    }
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselITest() {
            for (double nu = -16; nu <= 16; nu += 0.25) {
                Console.WriteLine(nu);

                foreach (double r in new double[] { 0, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64 }) {
                    foreach (double i in new double[] { 0, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64 }) {
                        if (r == 0 && i == 0) {
                            continue;
                        }

                        Complex<Pow2.N4> z = (r, i);

                        Console.WriteLine($"{z}: {BesselN4.BesselI(nu, z)}");
                    }
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselKTest() {
            for (double nu = 0; nu <= 16; nu += 0.25) {
                Console.WriteLine(nu);

                foreach (double r in new double[] { 0, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64 }) {
                    foreach (double i in new double[] { 0, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64 }) {
                        if (r == 0 && i == 0) {
                            continue;
                        }

                        Complex<Pow2.N4> z = (r, i);

                        Console.WriteLine($"{z}: {BesselN4.BesselK(nu, z)}");
                    }
                }

                Console.WriteLine(string.Empty);
            }
        }
    }
}