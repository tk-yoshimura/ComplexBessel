using ComplexBessel;
using DDoubleOptimizedBessel;
using DoubleDouble;
using MultiPrecision;
using MultiPrecisionComplex;

namespace DDoubleOptimizedBesselTests {
    [TestClass()]
    public class RealBesselExpectedTests {
        [TestMethod()]
        public void BesselJExpectedTest() {
            IEnumerable<string> filepaths = Directory.EnumerateFiles("../../../../mpmath_expected/besselj_real/", "*.csv");

            foreach (string filepath in filepaths) {
                string filename = filepath.Split("/").Last();
                string[] seps = filename.Split("_");
                ddouble x = ddouble.Parse(seps[2][1..^4]);
                using StreamReader sr = new(filepath);

                for (double nu = -16; nu <= 16; nu += 0.25) {
                    string? line = sr.ReadLine();

                    if (string.IsNullOrEmpty(line)) {
                        break;
                    }

                    ddouble expected = line;

                    ddouble actual = RealBessel.BesselJ(nu, x);
                    ddouble err = ddouble.Abs((expected - actual) / expected);

                    Console.WriteLine($"{nu}, {x}, {err:e10}");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    Assert.IsTrue(err < 4e-27, $"\n{nu}, {x}\n{expected}\n{actual}\n{err}");
                }
            }
        }

        [TestMethod()]
        public void BesselYExpectedTest() {
            IEnumerable<string> filepaths = Directory.EnumerateFiles("../../../../mpmath_expected/bessely_real/", "*.csv");

            foreach (string filepath in filepaths) {
                string filename = filepath.Split("/").Last();
                string[] seps = filename.Split("_");
                ddouble x = ddouble.Parse(seps[2][1..^4]);
                using StreamReader sr = new(filepath);

                for (double nu = -16; nu <= 16; nu += 0.25) {
                    string? line = sr.ReadLine();

                    if (string.IsNullOrEmpty(line)) {
                        break;
                    }

                    ddouble expected = line;

                    ddouble actual = RealBessel.BesselY(nu, x);
                    ddouble err = ddouble.Abs((expected - actual) / expected);

                    Console.WriteLine($"{nu}, {x}, {err:e10}");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    Assert.IsTrue(err < 4e-27, $"\n{nu}, {x}\n{expected}\n{actual}\n{err}");
                }
            }
        }

        [TestMethod()]
        public void BesselIExpectedTest() {
            IEnumerable<string> filepaths = Directory.EnumerateFiles("../../../../mpmath_expected/besseli_real/", "*.csv");

            foreach (string filepath in filepaths) {
                string filename = filepath.Split("/").Last();
                string[] seps = filename.Split("_");
                ddouble x = ddouble.Parse(seps[2][1..^4]);
                using StreamReader sr = new(filepath);

                for (double nu = -16; nu <= 16; nu += 0.25) {
                    string? line = sr.ReadLine();

                    if (string.IsNullOrEmpty(line)) {
                        break;
                    }

                    ddouble expected = line;

                    ddouble actual = RealBessel.BesselI(nu, x);
                    ddouble err = ddouble.Abs((expected - actual) / expected);

                    Console.WriteLine($"{nu}, {x}, {err:e10}");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    Assert.IsTrue(err < 2e-30, $"\n{nu}, {x}\n{expected}\n{actual}\n{err}");
                }
            }
        }

        [TestMethod()]
        public void BesselKExpectedTest() {
            IEnumerable<string> filepaths = Directory.EnumerateFiles("../../../../mpmath_expected/besselk_real/", "*.csv");

            foreach (string filepath in filepaths) {
                string filename = filepath.Split("/").Last();
                string[] seps = filename.Split("_");
                ddouble x = ddouble.Parse(seps[2][1..^4]);
                using StreamReader sr = new(filepath);

                for (double nu = -16; nu <= 16; nu += 0.25) {
                    string? line = sr.ReadLine();

                    if (string.IsNullOrEmpty(line)) {
                        break;
                    }

                    ddouble expected = line;

                    ddouble actual = RealBessel.BesselK(nu, x);
                    ddouble err = ddouble.Abs((expected - actual) / expected);

                    Console.WriteLine($"{nu}, {x}, {err:e10}");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    Assert.IsTrue(err < 8e-30, $"\n{nu}, {x}\n{expected}\n{actual}\n{err}");
                }
            }
        }
    }
}