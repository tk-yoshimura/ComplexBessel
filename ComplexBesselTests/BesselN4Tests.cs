using ComplexBessel;
using MultiPrecision;
using MultiPrecisionComplex;

namespace ComplexBesselTests {
    [TestClass()]
    public class BesselN4Tests {
        [TestMethod()]
        public void BesselJTest() {
            BesselN4.BesselJ(-15.5, (8.5, 7.5));

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

        [TestMethod()]
        public void BesselJExpectedTest() {
            IEnumerable<string> filepaths = Directory.EnumerateFiles("../../../../mpmath_expected/besselj/", "*.csv");

            foreach (string filepath in filepaths) {
                string filename = filepath.Split("/").Last();
                string[] seps = filename.Split("_");
                MultiPrecision<Pow2.N4> nu = MultiPrecision<Pow2.N4>.Parse(seps[1][2..^4]);
                using StreamReader sr = new(filepath);

                sr.ReadLine();

                while (!sr.EndOfStream) {
                    string? line = sr.ReadLine();

                    if (string.IsNullOrEmpty(line)) {
                        break;
                    }

                    string[] items = line.Split(",");

                    Complex<Pow2.N4> z = (items[0], items[1]);

                    if (!Complex<Pow2.N4>.TryParse(items[2], null, out Complex<Pow2.N4>? expected)) {
                        continue;
                    }

                    Complex<Pow2.N4> actual = BesselN4.BesselJ(nu, z);
                    MultiPrecision<Pow2.N4> err = (expected - actual).Magnitude / expected.Magnitude;

                    Console.WriteLine($"{nu}, {z}");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    Assert.IsTrue(err < 1e-34, $"\n{nu}, {z}\n{expected}\n{actual}\n{err}");
                }
            }
        }

        [TestMethod()]
        public void BesselYExpectedTest() {
            IEnumerable<string> filepaths = Directory.EnumerateFiles("../../../../mpmath_expected/bessely/", "*.csv");

            foreach (string filepath in filepaths) {
                string filename = filepath.Split("/").Last();
                string[] seps = filename.Split("_");
                MultiPrecision<Pow2.N4> nu = MultiPrecision<Pow2.N4>.Parse(seps[1][2..^4]);
                using StreamReader sr = new(filepath);

                sr.ReadLine();

                while (!sr.EndOfStream) {
                    string? line = sr.ReadLine();

                    if (string.IsNullOrEmpty(line)) {
                        break;
                    }

                    string[] items = line.Split(",");

                    Complex<Pow2.N4> z = (items[0], items[1]);

                    if (!Complex<Pow2.N4>.TryParse(items[2], null, out Complex<Pow2.N4>? expected)) {
                        continue;
                    }

                    Complex<Pow2.N4> actual = BesselN4.BesselY(nu, z);
                    MultiPrecision<Pow2.N4> err = (expected - actual).Magnitude / expected.Magnitude;

                    Console.WriteLine($"{nu}, {z}");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    if (z.R >= z.I * 2 || MultiPrecision<Pow2.N4>.Abs(nu) <= 2 || !MultiPrecision<Pow2.N4>.IsInteger(nu)) {
                        Assert.IsTrue(err < 1e-33, $"\n{nu}, {z}\n{expected}\n{actual}\n{err}");
                    }
                    else if (z.R >= z.I || MultiPrecision<Pow2.N4>.Abs(nu) <= 4) {
                        Assert.IsTrue(err < 4e-32, $"\n{nu}, {z}\n{expected}\n{actual}\n{err}");
                    }
                    else if (z.R >= z.I / 2 || MultiPrecision<Pow2.N4>.Abs(nu) <= 8) {
                        Assert.IsTrue(err < 2e-30, $"\n{nu}, {z}\n{expected}\n{actual}\n{err}");
                    }
                    else {
                        Assert.IsTrue(err < 8e-29, $"\n{nu}, {z}\n{expected}\n{actual}\n{err}");
                    }
                }
            }
        }

        [TestMethod()]
        public void BesselIExpectedTest() {
            IEnumerable<string> filepaths = Directory.EnumerateFiles("../../../../mpmath_expected/besseli/", "*.csv");

            foreach (string filepath in filepaths) {
                string filename = filepath.Split("/").Last();
                string[] seps = filename.Split("_");
                MultiPrecision<Pow2.N4> nu = MultiPrecision<Pow2.N4>.Parse(seps[1][2..^4]);
                using StreamReader sr = new(filepath);

                sr.ReadLine();

                while (!sr.EndOfStream) {
                    string? line = sr.ReadLine();

                    if (string.IsNullOrEmpty(line)) {
                        break;
                    }

                    string[] items = line.Split(",");

                    Complex<Pow2.N4> z = (items[0], items[1]);

                    if (!Complex<Pow2.N4>.TryParse(items[2], null, out Complex<Pow2.N4>? expected)) {
                        continue;
                    }

                    Complex<Pow2.N4> actual = BesselN4.BesselI(nu, z);
                    MultiPrecision<Pow2.N4> err = (expected - actual).Magnitude / expected.Magnitude;

                    Console.WriteLine($"{nu}, {z}");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    Assert.IsTrue(err < 1e-34, $"\n{nu}, {z}\n{expected}\n{actual}\n{err}");
                }
            }
        }
    }
}