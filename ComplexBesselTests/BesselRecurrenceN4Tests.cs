using ComplexBessel;
using MultiPrecision;
using MultiPrecisionComplex;

namespace ComplexBesselTests {
    [TestClass()]
    public class BesselRecurrenceN4Tests {
        [TestMethod()]
        public void BesselJRecurrenceExpectedTest() {
            IEnumerable<string> filepaths = Directory.EnumerateFiles("../../../../mpmath_expected/besselj_complex_nuplus16/", "*.csv");

            foreach (string filepath in filepaths) {
                string filename = filepath.Split("/").Last();
                string[] seps = filename.Split("_");
                MultiPrecision<Pow2.N4> nu = MultiPrecision<Pow2.N4>.Parse(seps[3][2..^4]);
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
        public void BesselYRecurrenceExpectedTest() {
            IEnumerable<string> filepaths = Directory.EnumerateFiles("../../../../mpmath_expected/bessely_complex_nuplus16/", "*.csv");

            foreach (string filepath in filepaths) {
                string filename = filepath.Split("/").Last();
                string[] seps = filename.Split("_");
                MultiPrecision<Pow2.N4> nu = MultiPrecision<Pow2.N4>.Parse(seps[3][2..^4]);
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

                    Assert.IsTrue(err < 1e-34, $"\n{nu}, {z}\n{expected}\n{actual}\n{err}");
                }
            }
        }

        [TestMethod()]
        public void BesselIRecurrenceExpectedTest() {
            IEnumerable<string> filepaths = Directory.EnumerateFiles("../../../../mpmath_expected/besseli_complex_nuplus16/", "*.csv");

            foreach (string filepath in filepaths) {
                string filename = filepath.Split("/").Last();
                string[] seps = filename.Split("_");
                MultiPrecision<Pow2.N4> nu = MultiPrecision<Pow2.N4>.Parse(seps[3][2..^4]);
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

        [TestMethod()]
        public void BesselKRecurrenceExpectedTest() {
            IEnumerable<string> filepaths = Directory.EnumerateFiles("../../../../mpmath_expected/besselk_complex_nuplus16/", "*.csv");

            foreach (string filepath in filepaths) {
                string filename = filepath.Split("/").Last();
                string[] seps = filename.Split("_");
                MultiPrecision<Pow2.N4> nu = MultiPrecision<Pow2.N4>.Parse(seps[3][2..^4]);
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

                    Complex<Pow2.N4> actual = BesselN4.BesselK(nu, z);
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