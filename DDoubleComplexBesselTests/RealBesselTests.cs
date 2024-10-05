using ComplexBessel;
using DDoubleComplexBessel;
using DoubleDouble;
using DoubleDoubleComplex;

namespace DDoubleComplexBesselTests {
    [TestClass()]
    public class RealBesselTests {
        [TestMethod()]
        public void BesselJTest() {
            for (double nu = -16; nu <= 16; nu += 0.25) {
                Console.WriteLine(nu);

                for (double r = 0.5; r <= 42; r += 0.5) {
                    Complex expected = BesselN4.BesselJ(nu, r).ToString();
                    Complex actual = Bessel.BesselJ(nu, r);

                    ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                    Console.WriteLine($"{nu}, {r}, {err:e4}");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    Assert.IsTrue(err < 2e-27, $"\n{nu}, {r}\n{expected}\n{actual}\n{err:e4}");
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselYTest() {
            for (double nu = -16; nu <= 16; nu += 0.25) {
                Console.WriteLine(nu);

                for (double r = 0.5; r <= 42; r += 0.5) {
                    Complex expected = BesselN4.BesselY(nu, r).ToString();
                    Complex actual = Bessel.BesselY(nu, r);

                    ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                    Console.WriteLine($"{nu}, {r}, {err:e4}");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    Assert.IsTrue(err < 2e-27, $"\n{nu}, {r}\n{expected}\n{actual}\n{err:e4}");
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselITest() {
            for (double nu = -16; nu <= 16; nu += 0.25) {
                Console.WriteLine(nu);

                for (double r = 0.125; r <= 42; r += 0.125) {
                    Complex expected = BesselN4.BesselI(nu, r).ToString();
                    Complex actual = Bessel.BesselI(nu, r);

                    ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                    Console.WriteLine($"{nu}, {r}, {err:e4}");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    Assert.IsTrue(err < 2e-27, $"\n{nu}, {r}\n{expected}\n{actual}\n{err:e4}");
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselKTest() {
            for (double nu = 0; nu <= 16; nu += 0.25) {
                Console.WriteLine(nu);

                for (double r = 0.125; r <= 42; r += 0.125) {
                    Complex expected = BesselN4.BesselK(nu, r).ToString();
                    Complex actual = Bessel.BesselK(nu, r);

                    ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                    Console.WriteLine($"{nu}, {r}, {err:e4}");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    Assert.IsTrue(err < 2e-27, $"\n{nu}, {r}\n{expected}\n{actual}\n{err:e4}");
                }

                Console.WriteLine(string.Empty);
            }
        }
    }
}