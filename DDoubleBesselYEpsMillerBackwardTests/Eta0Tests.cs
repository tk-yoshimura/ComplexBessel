using DoubleDouble;
using MultiPrecision;

namespace DDoubleBesselYEpsMillerBackwardTests {
    [TestClass]
    public class Eta0Tests {

        [TestMethod]
        public void AlphaPlusExpM10Test() {
            double eps = double.ScaleB(1, -10);

            for (double x = 1; x <= 64; x += 0.5) {
                ddouble eta0_eps = DDoubleCoef.BesselYEta0Xi1Eps(eps, x).eta0;
                ddouble eta0_raw = DDoubleCoef.BesselYEta0Xi1(eps, x).eta0;
                ddouble eta0_mp = MPCoef<Pow2.N16>.BesselYEta0Xi1(eps, x).eta0.ToString();

                ddouble err_eps = ddouble.Abs((eta0_eps - eta0_mp) / eta0_mp);
                ddouble err_raw = ddouble.Abs((eta0_raw - eta0_mp) / eta0_mp);

                Console.WriteLine(x);

                Console.WriteLine($"exp: {eta0_mp}");
                Console.WriteLine($"eps: {eta0_eps}");
                Console.WriteLine($"raw: {eta0_raw}");

                Console.WriteLine($"eps: {err_eps:e4}");
                Console.WriteLine($"raw: {err_raw:e4}");

                Console.WriteLine("");
            }
        }

        [TestMethod]
        public void AlphaPlusExpM12Test() {
            double eps = double.ScaleB(1, -12);

            for (double x = 1; x <= 64; x += 0.5) {
                ddouble eta0_eps = DDoubleCoef.BesselYEta0Xi1Eps(eps, x).eta0;
                ddouble eta0_raw = DDoubleCoef.BesselYEta0Xi1(eps, x).eta0;
                ddouble eta0_mp = MPCoef<Pow2.N16>.BesselYEta0Xi1(eps, x).eta0.ToString();

                ddouble err_eps = ddouble.Abs((eta0_eps - eta0_mp) / eta0_mp);
                ddouble err_raw = ddouble.Abs((eta0_raw - eta0_mp) / eta0_mp);

                Console.WriteLine(x);

                Console.WriteLine($"exp: {eta0_mp}");
                Console.WriteLine($"eps: {eta0_eps}");
                Console.WriteLine($"raw: {eta0_raw}");

                Console.WriteLine($"eps: {err_eps:e4}");
                Console.WriteLine($"raw: {err_raw:e4}");

                Console.WriteLine("");
            }
        }

        [TestMethod]
        public void AlphaPlusExpM16Test() {
            double eps = double.ScaleB(1, -16);

            for (double x = 1; x <= 64; x += 0.5) {
                ddouble eta0_eps = DDoubleCoef.BesselYEta0Xi1Eps(eps, x).eta0;
                ddouble eta0_raw = DDoubleCoef.BesselYEta0Xi1(eps, x).eta0;
                ddouble eta0_mp = MPCoef<Pow2.N16>.BesselYEta0Xi1(eps, x).eta0.ToString();

                ddouble err_eps = ddouble.Abs((eta0_eps - eta0_mp) / eta0_mp);
                ddouble err_raw = ddouble.Abs((eta0_raw - eta0_mp) / eta0_mp);

                Console.WriteLine(x);

                Console.WriteLine($"exp: {eta0_mp}");
                Console.WriteLine($"eps: {eta0_eps}");
                Console.WriteLine($"raw: {eta0_raw}");

                Console.WriteLine($"eps: {err_eps:e4}");
                Console.WriteLine($"raw: {err_raw:e4}");

                Console.WriteLine("");
            }
        }

        [TestMethod]
        public void AlphaPlusExpM20Test() {
            double eps = double.ScaleB(1, -20);

            for (double x = 1; x <= 64; x += 0.5) {
                ddouble eta0_eps = DDoubleCoef.BesselYEta0Xi1Eps(eps, x).eta0;
                ddouble eta0_raw = DDoubleCoef.BesselYEta0Xi1(eps, x).eta0;
                ddouble eta0_mp = MPCoef<Pow2.N16>.BesselYEta0Xi1(eps, x).eta0.ToString();

                ddouble err_eps = ddouble.Abs((eta0_eps - eta0_mp) / eta0_mp);
                ddouble err_raw = ddouble.Abs((eta0_raw - eta0_mp) / eta0_mp);

                Console.WriteLine(x);

                Console.WriteLine($"exp: {eta0_mp}");
                Console.WriteLine($"eps: {eta0_eps}");
                Console.WriteLine($"raw: {eta0_raw}");

                Console.WriteLine($"eps: {err_eps:e4}");
                Console.WriteLine($"raw: {err_raw:e4}");

                Console.WriteLine("");
            }
        }

        [TestMethod]
        public void AlphaPlusExpM24Test() {
            double eps = double.ScaleB(1, -24);

            for (double x = 1; x <= 64; x += 0.5) {
                ddouble eta0_eps = DDoubleCoef.BesselYEta0Xi1Eps(eps, x).eta0;
                ddouble eta0_raw = DDoubleCoef.BesselYEta0Xi1(eps, x).eta0;
                ddouble eta0_mp = MPCoef<Pow2.N16>.BesselYEta0Xi1(eps, x).eta0.ToString();

                ddouble err_eps = ddouble.Abs((eta0_eps - eta0_mp) / eta0_mp);
                ddouble err_raw = ddouble.Abs((eta0_raw - eta0_mp) / eta0_mp);

                Console.WriteLine(x);

                Console.WriteLine($"exp: {eta0_mp}");
                Console.WriteLine($"eps: {eta0_eps}");
                Console.WriteLine($"raw: {eta0_raw}");

                Console.WriteLine($"eps: {err_eps:e4}");
                Console.WriteLine($"raw: {err_raw:e4}");

                Console.WriteLine("");
            }
        }

        [TestMethod]
        public void AlphaPlusExpM28Test() {
            double eps = double.ScaleB(1, -28);

            for (double x = 1; x <= 64; x += 0.5) {
                ddouble eta0_eps = DDoubleCoef.BesselYEta0Xi1Eps(eps, x).eta0;
                ddouble eta0_raw = DDoubleCoef.BesselYEta0Xi1(eps, x).eta0;
                ddouble eta0_mp = MPCoef<Pow2.N16>.BesselYEta0Xi1(eps, x).eta0.ToString();

                ddouble err_eps = ddouble.Abs((eta0_eps - eta0_mp) / eta0_mp);
                ddouble err_raw = ddouble.Abs((eta0_raw - eta0_mp) / eta0_mp);

                Console.WriteLine(x);

                Console.WriteLine($"exp: {eta0_mp}");
                Console.WriteLine($"eps: {eta0_eps}");
                Console.WriteLine($"raw: {eta0_raw}");

                Console.WriteLine($"eps: {err_eps:e4}");
                Console.WriteLine($"raw: {err_raw:e4}");

                Console.WriteLine("");
            }
        }

        [TestMethod]
        public void AlphaPlusExpM32Test() {
            double eps = double.ScaleB(1, -32);

            for (double x = 1; x <= 64; x += 0.5) {
                ddouble eta0_eps = DDoubleCoef.BesselYEta0Xi1Eps(eps, x).eta0;
                ddouble eta0_raw = DDoubleCoef.BesselYEta0Xi1(eps, x).eta0;
                ddouble eta0_mp = MPCoef<Pow2.N16>.BesselYEta0Xi1(eps, x).eta0.ToString();

                ddouble err_eps = ddouble.Abs((eta0_eps - eta0_mp) / eta0_mp);
                ddouble err_raw = ddouble.Abs((eta0_raw - eta0_mp) / eta0_mp);

                Console.WriteLine(x);

                Console.WriteLine($"exp: {eta0_mp}");
                Console.WriteLine($"eps: {eta0_eps}");
                Console.WriteLine($"raw: {eta0_raw}");

                Console.WriteLine($"eps: {err_eps:e4}");
                Console.WriteLine($"raw: {err_raw:e4}");

                Console.WriteLine("");
            }
        }

        [TestMethod]
        public void AlphaMinusExpM10Test() {
            double eps = -double.ScaleB(1, -10);

            for (double x = 1; x <= 64; x += 0.5) {
                ddouble eta0_eps = DDoubleCoef.BesselYEta0Xi1Eps(eps, x).eta0;
                ddouble eta0_raw = DDoubleCoef.BesselYEta0Xi1(eps, x).eta0;
                ddouble eta0_mp = MPCoef<Pow2.N16>.BesselYEta0Xi1(eps, x).eta0.ToString();

                ddouble err_eps = ddouble.Abs((eta0_eps - eta0_mp) / eta0_mp);
                ddouble err_raw = ddouble.Abs((eta0_raw - eta0_mp) / eta0_mp);

                Console.WriteLine(x);

                Console.WriteLine($"exp: {eta0_mp}");
                Console.WriteLine($"eps: {eta0_eps}");
                Console.WriteLine($"raw: {eta0_raw}");

                Console.WriteLine($"eps: {err_eps:e4}");
                Console.WriteLine($"raw: {err_raw:e4}");

                Console.WriteLine("");
            }
        }

        [TestMethod]
        public void AlphaMinusExpM12Test() {
            double eps = -double.ScaleB(1, -12);

            for (double x = 1; x <= 64; x += 0.5) {
                ddouble eta0_eps = DDoubleCoef.BesselYEta0Xi1Eps(eps, x).eta0;
                ddouble eta0_raw = DDoubleCoef.BesselYEta0Xi1(eps, x).eta0;
                ddouble eta0_mp = MPCoef<Pow2.N16>.BesselYEta0Xi1(eps, x).eta0.ToString();

                ddouble err_eps = ddouble.Abs((eta0_eps - eta0_mp) / eta0_mp);
                ddouble err_raw = ddouble.Abs((eta0_raw - eta0_mp) / eta0_mp);

                Console.WriteLine(x);

                Console.WriteLine($"exp: {eta0_mp}");
                Console.WriteLine($"eps: {eta0_eps}");
                Console.WriteLine($"raw: {eta0_raw}");

                Console.WriteLine($"eps: {err_eps:e4}");
                Console.WriteLine($"raw: {err_raw:e4}");

                Console.WriteLine("");
            }
        }

        [TestMethod]
        public void AlphaMinusExpM16Test() {
            double eps = -double.ScaleB(1, -16);

            for (double x = 1; x <= 64; x += 0.5) {
                ddouble eta0_eps = DDoubleCoef.BesselYEta0Xi1Eps(eps, x).eta0;
                ddouble eta0_raw = DDoubleCoef.BesselYEta0Xi1(eps, x).eta0;
                ddouble eta0_mp = MPCoef<Pow2.N16>.BesselYEta0Xi1(eps, x).eta0.ToString();

                ddouble err_eps = ddouble.Abs((eta0_eps - eta0_mp) / eta0_mp);
                ddouble err_raw = ddouble.Abs((eta0_raw - eta0_mp) / eta0_mp);

                Console.WriteLine(x);

                Console.WriteLine($"exp: {eta0_mp}");
                Console.WriteLine($"eps: {eta0_eps}");
                Console.WriteLine($"raw: {eta0_raw}");

                Console.WriteLine($"eps: {err_eps:e4}");
                Console.WriteLine($"raw: {err_raw:e4}");

                Console.WriteLine("");
            }
        }

        [TestMethod]
        public void AlphaMinusExpM20Test() {
            double eps = -double.ScaleB(1, -20);

            for (double x = 1; x <= 64; x += 0.5) {
                ddouble eta0_eps = DDoubleCoef.BesselYEta0Xi1Eps(eps, x).eta0;
                ddouble eta0_raw = DDoubleCoef.BesselYEta0Xi1(eps, x).eta0;
                ddouble eta0_mp = MPCoef<Pow2.N16>.BesselYEta0Xi1(eps, x).eta0.ToString();

                ddouble err_eps = ddouble.Abs((eta0_eps - eta0_mp) / eta0_mp);
                ddouble err_raw = ddouble.Abs((eta0_raw - eta0_mp) / eta0_mp);

                Console.WriteLine(x);

                Console.WriteLine($"exp: {eta0_mp}");
                Console.WriteLine($"eps: {eta0_eps}");
                Console.WriteLine($"raw: {eta0_raw}");

                Console.WriteLine($"eps: {err_eps:e4}");
                Console.WriteLine($"raw: {err_raw:e4}");

                Console.WriteLine("");
            }
        }

        [TestMethod]
        public void AlphaMinusExpM24Test() {
            double eps = -double.ScaleB(1, -24);

            for (double x = 1; x <= 64; x += 0.5) {
                ddouble eta0_eps = DDoubleCoef.BesselYEta0Xi1Eps(eps, x).eta0;
                ddouble eta0_raw = DDoubleCoef.BesselYEta0Xi1(eps, x).eta0;
                ddouble eta0_mp = MPCoef<Pow2.N16>.BesselYEta0Xi1(eps, x).eta0.ToString();

                ddouble err_eps = ddouble.Abs((eta0_eps - eta0_mp) / eta0_mp);
                ddouble err_raw = ddouble.Abs((eta0_raw - eta0_mp) / eta0_mp);

                Console.WriteLine(x);

                Console.WriteLine($"exp: {eta0_mp}");
                Console.WriteLine($"eps: {eta0_eps}");
                Console.WriteLine($"raw: {eta0_raw}");

                Console.WriteLine($"eps: {err_eps:e4}");
                Console.WriteLine($"raw: {err_raw:e4}");

                Console.WriteLine("");
            }
        }

        [TestMethod]
        public void AlphaMinusExpM28Test() {
            double eps = -double.ScaleB(1, -28);

            for (double x = 1; x <= 64; x += 0.5) {
                ddouble eta0_eps = DDoubleCoef.BesselYEta0Xi1Eps(eps, x).eta0;
                ddouble eta0_raw = DDoubleCoef.BesselYEta0Xi1(eps, x).eta0;
                ddouble eta0_mp = MPCoef<Pow2.N16>.BesselYEta0Xi1(eps, x).eta0.ToString();

                ddouble err_eps = ddouble.Abs((eta0_eps - eta0_mp) / eta0_mp);
                ddouble err_raw = ddouble.Abs((eta0_raw - eta0_mp) / eta0_mp);

                Console.WriteLine(x);

                Console.WriteLine($"exp: {eta0_mp}");
                Console.WriteLine($"eps: {eta0_eps}");
                Console.WriteLine($"raw: {eta0_raw}");

                Console.WriteLine($"eps: {err_eps:e4}");
                Console.WriteLine($"raw: {err_raw:e4}");

                Console.WriteLine("");
            }
        }

        [TestMethod]
        public void AlphaMinusExpM32Test() {
            double eps = -double.ScaleB(1, -32);

            for (double x = 1; x <= 64; x += 0.5) {
                ddouble eta0_eps = DDoubleCoef.BesselYEta0Xi1Eps(eps, x).eta0;
                ddouble eta0_raw = DDoubleCoef.BesselYEta0Xi1(eps, x).eta0;
                ddouble eta0_mp = MPCoef<Pow2.N16>.BesselYEta0Xi1(eps, x).eta0.ToString();

                ddouble err_eps = ddouble.Abs((eta0_eps - eta0_mp) / eta0_mp);
                ddouble err_raw = ddouble.Abs((eta0_raw - eta0_mp) / eta0_mp);

                Console.WriteLine(x);

                Console.WriteLine($"exp: {eta0_mp}");
                Console.WriteLine($"eps: {eta0_eps}");
                Console.WriteLine($"raw: {eta0_raw}");

                Console.WriteLine($"eps: {err_eps:e4}");
                Console.WriteLine($"raw: {err_raw:e4}");

                Console.WriteLine("");
            }
        }

        [TestMethod]
        public void AccuracyTest() {

            for (int exp = -32; exp <= -10; exp++) {
                double eps = -double.ScaleB(1, exp);

                Console.WriteLine($"eps: 2^{eps}");

                for (double x = 1; x <= 64; x += 0.5) {
                    ddouble eta0_eps = DDoubleCoef.BesselYEta0Xi1Eps(eps, x).eta0;
                    ddouble eta0_raw = DDoubleCoef.BesselYEta0Xi1(eps, x).eta0;
                    ddouble eta0_mp = MPCoef<Pow2.N16>.BesselYEta0Xi1(eps, x).eta0.ToString();

                    ddouble err_eps = ddouble.Abs((eta0_eps - eta0_mp) / eta0_mp);
                    ddouble err_raw = ddouble.Abs((eta0_raw - eta0_mp) / eta0_mp);

                    Console.WriteLine(x);

                    Console.WriteLine($"exp: {eta0_mp}");
                    Console.WriteLine($"eps: {eta0_eps}");
                    Console.WriteLine($"raw: {eta0_raw}");

                    Console.WriteLine($"eps: {err_eps:e4}");
                    Console.WriteLine($"raw: {err_raw:e4}");

                    Assert.IsTrue(err_eps < 1e-27 || err_raw < 1e-27, $"2^{exp}, {x}, {err_eps: e4}, {err_raw: e4}");

                    Console.WriteLine("");
                }

                Console.WriteLine($"eps: -2^{eps}");

                for (double x = 1; x <= 64; x += 0.5) {
                    ddouble eta0_eps = DDoubleCoef.BesselYEta0Xi1Eps(-eps, x).eta0;
                    ddouble eta0_raw = DDoubleCoef.BesselYEta0Xi1(-eps, x).eta0;
                    ddouble eta0_mp = MPCoef<Pow2.N16>.BesselYEta0Xi1(-eps, x).eta0.ToString();

                    ddouble err_eps = ddouble.Abs((eta0_eps - eta0_mp) / eta0_mp);
                    ddouble err_raw = ddouble.Abs((eta0_raw - eta0_mp) / eta0_mp);

                    Console.WriteLine(x);

                    Console.WriteLine($"exp: {eta0_mp}");
                    Console.WriteLine($"eps: {eta0_eps}");
                    Console.WriteLine($"raw: {eta0_raw}");

                    Console.WriteLine($"eps: {err_eps:e4}");
                    Console.WriteLine($"raw: {err_raw:e4}");

                    Assert.IsTrue(err_eps < 1e-26 || err_raw < 1e-26, $"-2^{exp}, {x}, {err_eps: e4}, {err_raw: e4}");

                    Console.WriteLine("");
                }
            }
        }

        [TestMethod]
        public void R2Test() {
            for (int exp = -32; exp <= -10; exp++) {
                double eps = -double.ScaleB(1, exp);

                Console.WriteLine($"eps: 2^{eps}");

                for (double x = 1; x <= 64; x += 0.5) {
                    ddouble eta0_eps = DDoubleCoef.BesselYEta0Xi1Eps(eps, x).eta0;
                    ddouble eta0_r2 = DDoubleCoef.BesselYEta0Xi1EpsR2(eps, x).eta0;

                    Assert.AreEqual(eta0_eps, eta0_r2);
                }

                Console.WriteLine($"eps: -2^{eps}");

                for (double x = 1; x <= 64; x += 0.5) {
                    ddouble eta0_eps = DDoubleCoef.BesselYEta0Xi1Eps(eps, x).eta0;
                    ddouble eta0_r2 = DDoubleCoef.BesselYEta0Xi1EpsR2(eps, x).eta0;

                    Assert.AreEqual(eta0_eps, eta0_r2);
                }
            }
        }
    }
}