using ComplexBessel;
using DDoubleComplexBessel;
using DoubleDouble;
using DoubleDoubleComplex;

namespace DDoubleComplexBesselBenchmark {
    [TestClass]
    public class MillerBackwardIteration {
        [TestMethod()]
        public void BesselJTest() {
            for (double nu = -16; nu <= 16; nu += 0.25) {
                using StreamWriter sw = new($"../../../../results_ddouble/besselj_nu{nu}_millerbackward_iteration.csv");
                sw.WriteLine("r,i,m");

                for (double r = 5; r <= 40; r += 0.25) {

                    for (double i = 0; i <= 6; i += 0.125) {
                        if (r * r + i * i > 40 * 40) {
                            continue;
                        }

                        Complex expected = BesselN4.BesselJ(nu, (r, i)).ToString();

                        int m;
                        for (m = 20; m <= 256; m += 16) {
                            Complex actual = MillerBackward.BesselJKernel(nu, (r, i), m);

                            ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                            if (err < 3.16e-30) {
                                break;
                            }
                        }

                        for (m = int.Max(20, m - 16); m <= 256; m += 2) {
                            Complex actual = MillerBackward.BesselJKernel(nu, (r, i), m);

                            ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                            if (err < 3.16e-30) {
                                sw.WriteLine($"{r},{i},{m}");
                                break;
                            }
                        }
                    }
                }
            }
        }

        [TestMethod()]
        public void BesselYTest() {
            for (double nu = -16; nu <= 16; nu += 0.25) {
                using StreamWriter sw = new($"../../../../results_ddouble/bessely_nu{nu}_millerbackward_iteration.csv");
                sw.WriteLine("r,i,m");

                for (double r = 5; r <= 40; r += 0.25) {

                    for (double i = 0; i <= 6; i += 0.125) {
                        if (r * r + i * i > 40 * 40) {
                            continue;
                        }

                        Complex expected = BesselN4.BesselY(nu, (r, i)).ToString();

                        int m;
                        for (m = 20; m <= 256; m += 16) {
                            Complex actual = MillerBackward.BesselYKernel(nu, (r, i), m);

                            ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                            if (err < 3.16e-30) {
                                break;
                            }
                        }

                        for (m = int.Max(20, m - 16); m <= 256; m += 2) {
                            Complex actual = MillerBackward.BesselYKernel(nu, (r, i), m);

                            ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                            if (err < 3.16e-30) {
                                sw.WriteLine($"{r},{i},{m}");
                                break;
                            }
                        }
                    }
                }
            }
        }

        [TestMethod()]
        public void BesselITest() {
            for (double nu = -16; nu <= 16; nu += 0.25) {
                using StreamWriter sw = new($"../../../../results_ddouble/besseli_nu{nu}_millerbackward_iteration.csv");
                sw.WriteLine("r,i,m");

                for (double r = 0; r <= 40; r += 0.25) {
                    for (double i = 6; i <= 40; i += 0.25) {

                        if (r * r + i * i > 40 * 40) {
                            continue;
                        }

                        Complex expected = BesselN4.BesselI(nu, (r, i)).ToString();

                        int m;
                        for (m = 20; m <= 256; m += 16) {
                            Complex actual = MillerBackward.BesselIKernel(nu, (r, i), m);

                            ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                            if (err < 3.16e-30) {
                                break;
                            }
                        }

                        for (m = int.Max(20, m - 16); m <= 256; m += 2) {
                            Complex actual = MillerBackward.BesselIKernel(nu, (r, i), m);

                            ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                            if (err < 3.16e-30) {
                                sw.WriteLine($"{r},{i},{m}");
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
}