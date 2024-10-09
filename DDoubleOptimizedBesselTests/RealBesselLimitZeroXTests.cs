using DoubleDouble;

namespace DDoubleOptimizedBesselTests {
    [TestClass()]
    public class RealBesselLimitZeroXTests {
        [TestMethod]
        public void BesselJLimitZeroXTest() {
            ddouble eps = ddouble.Ldexp(1, -4);

            for (ddouble nu = -16; nu <= 16; nu += 1d / 32) {
                Console.WriteLine($"nu={nu}");

                ddouble y0 = DDoubleOptimizedBessel.RealBessel.BesselJ(nu, eps * 4);
                ddouble y1 = DDoubleOptimizedBessel.RealBessel.BesselJ(nu, eps * 2);

                if (ddouble.IsInfinity(y0) || ddouble.IsInfinity(y1)) {
                    continue;
                }

                for (ddouble x = eps; x > 0; x /= 2) {
                    ddouble y2 = DDoubleOptimizedBessel.RealBessel.BesselJ(nu, x);

                    Console.WriteLine($"x = {x}\ny = {y2}");

                    Assert.IsFalse(ddouble.IsNaN(y2), $"{nu},{x}");

                    if (y1 > 1e-300) {
                        if (ddouble.Abs(y0) >= ddouble.Abs(y1)) {
                            Assert.IsTrue(ddouble.Abs(y1) >= ddouble.Abs(y2), $"{nu},{x}");
                        }
                        if (ddouble.Abs(y0) <= ddouble.Abs(y1)) {
                            Assert.IsTrue(ddouble.Abs(y1) <= ddouble.Abs(y2), $"{nu},{x}");
                        }
                    }

                    if (y2 == 0) {
                        break;
                    }

                    Assert.AreEqual(ddouble.Sign(y1), ddouble.Sign(y2), $"{nu},{x}");

                    if (ddouble.IsInfinity(y2)) {
                        break;
                    }

                    (y0, y1) = (y1, y2);
                }
            }

            for (ddouble nu = -16; nu <= 16; nu += 1d / 256) {
                for (ddouble x = ddouble.Ldexp(1, -950); x > 0; x /= 2) {
                    ddouble y = DDoubleOptimizedBessel.RealBessel.BesselJ(nu, x);

                    Assert.IsFalse(ddouble.IsNaN(y), $"{nu},{x}");
                }

                ddouble y_zero = DDoubleOptimizedBessel.RealBessel.BesselJ(nu, 0);

                Assert.IsFalse(ddouble.IsNaN(y_zero), $"{nu},0");
            }

            for (int n = -15; n <= 15; n++) {
                foreach (ddouble alpha in new[] {
                    -double.ScaleB(1, -12), -double.ScaleB(1, -24), -double.ScaleB(1, -48), -double.ScaleB(1, -96),
                    double.ScaleB(1, -12), double.ScaleB(1, -24), double.ScaleB(1, -48), double.ScaleB(1, -96), }) {

                    ddouble nu = n + alpha;

                    for (ddouble x = ddouble.Ldexp(1, -950); x > 0; x /= 2) {
                        ddouble y = DDoubleOptimizedBessel.RealBessel.BesselJ(nu, x);

                        Assert.IsFalse(ddouble.IsNaN(y), $"{nu},{x}");
                    }

                    ddouble y_zero = DDoubleOptimizedBessel.RealBessel.BesselJ(nu, 0);

                    Assert.IsFalse(ddouble.IsNaN(y_zero), $"{nu},0");
                }
            }
        }

        [TestMethod]
        public void BesselYLimitZeroXTest() {
            ddouble eps = ddouble.Ldexp(1, -4);

            for (ddouble nu = -16; nu <= 16; nu += 1d / 32) {
                Console.WriteLine($"nu={nu}");

                ddouble y0 = DDoubleOptimizedBessel.RealBessel.BesselY(nu, eps * 4);
                ddouble y1 = DDoubleOptimizedBessel.RealBessel.BesselY(nu, eps * 2);

                if (ddouble.IsInfinity(y0) || ddouble.IsInfinity(y1)) {
                    continue;
                }

                for (ddouble x = eps; x > 0; x /= 2) {
                    ddouble y2 = DDoubleOptimizedBessel.RealBessel.BesselY(nu, x);

                    Console.WriteLine($"x = {x}\ny = {y2}");

                    Assert.IsFalse(ddouble.IsNaN(y2), $"{nu},{x}");

                    if (y1 > 1e-300 && nu != -0.53125) { // ignore: nu = -0.53125
                        if (ddouble.Abs(y0) >= ddouble.Abs(y1)) {
                            Assert.IsTrue(ddouble.Abs(y1) >= ddouble.Abs(y2), $"{nu},{x}");
                        }
                        if (ddouble.Abs(y0) <= ddouble.Abs(y1)) {
                            Assert.IsTrue(ddouble.Abs(y1) <= ddouble.Abs(y2), $"{nu},{x}");
                        }
                    }

                    if (y2 == 0) {
                        break;
                    }

                    if (nu != -0.46875) { // ignore: nu = -0.46875
                        Assert.AreEqual(ddouble.Sign(y1), ddouble.Sign(y2), $"{nu},{x}");
                    }

                    if (ddouble.IsInfinity(y2)) {
                        break;
                    }

                    (y0, y1) = (y1, y2);
                }
            }

            for (ddouble nu = -16; nu <= 16; nu += 1d / 256) {
                for (ddouble x = ddouble.Ldexp(1, -950); x > 0; x /= 2) {
                    ddouble y = DDoubleOptimizedBessel.RealBessel.BesselY(nu, x);

                    Assert.IsFalse(ddouble.IsNaN(y), $"{nu},{x}");
                }

                ddouble y_zero = DDoubleOptimizedBessel.RealBessel.BesselY(nu, 0);

                Assert.IsFalse(ddouble.IsNaN(y_zero), $"{nu},0");
            }

            for (int n = -15; n <= 15; n++) {
                foreach (ddouble alpha in new[] {
                    -double.ScaleB(1, -12), -double.ScaleB(1, -24), -double.ScaleB(1, -48), -double.ScaleB(1, -96),
                    double.ScaleB(1, -12), double.ScaleB(1, -24), double.ScaleB(1, -48), double.ScaleB(1, -96), }) {

                    ddouble nu = n + alpha;

                    for (ddouble x = ddouble.Ldexp(1, -950); x > 0; x /= 2) {
                        ddouble y = DDoubleOptimizedBessel.RealBessel.BesselY(nu, x);

                        Assert.IsFalse(ddouble.IsNaN(y), $"{nu},{x}");
                    }

                    ddouble y_zero = DDoubleOptimizedBessel.RealBessel.BesselY(nu, 0);

                    Assert.IsFalse(ddouble.IsNaN(y_zero), $"{nu},0");
                }
            }
        }

        [TestMethod]
        public void BesselILimitZeroXTest() {
            ddouble eps = ddouble.Ldexp(1, -4);

            for (ddouble nu = -16; nu <= 16; nu += 1d / 32) {
                Console.WriteLine($"nu={nu}");

                ddouble y0 = DDoubleOptimizedBessel.RealBessel.BesselI(nu, eps * 4);
                ddouble y1 = DDoubleOptimizedBessel.RealBessel.BesselI(nu, eps * 2);

                if (ddouble.IsInfinity(y0) || ddouble.IsInfinity(y1)) {
                    continue;
                }

                for (ddouble x = eps; x > 0; x /= 2) {
                    ddouble y2 = DDoubleOptimizedBessel.RealBessel.BesselI(nu, x);

                    Console.WriteLine($"x = {x}\ny = {y2}");

                    Assert.IsFalse(ddouble.IsNaN(y2), $"{nu},{x}");

                    if (y1 > 1e-300) {
                        if (ddouble.Abs(y0) >= ddouble.Abs(y1)) {
                            Assert.IsTrue(ddouble.Abs(y1) >= ddouble.Abs(y2), $"{nu},{x}");
                        }
                        if (ddouble.Abs(y0) <= ddouble.Abs(y1)) {
                            Assert.IsTrue(ddouble.Abs(y1) <= ddouble.Abs(y2), $"{nu},{x}");
                        }
                    }

                    if (y2 == 0) {
                        break;
                    }

                    Assert.AreEqual(ddouble.Sign(y1), ddouble.Sign(y2), $"{nu},{x}");

                    if (ddouble.IsInfinity(y2)) {
                        break;
                    }

                    (y0, y1) = (y1, y2);
                }
            }

            for (ddouble nu = -16; nu <= 16; nu += 1d / 256) {
                for (ddouble x = ddouble.Ldexp(1, -950); x > 0; x /= 2) {
                    ddouble y = DDoubleOptimizedBessel.RealBessel.BesselI(nu, x);

                    Assert.IsFalse(ddouble.IsNaN(y), $"{nu},{x}");
                }

                ddouble y_zero = DDoubleOptimizedBessel.RealBessel.BesselI(nu, 0);

                Assert.IsFalse(ddouble.IsNaN(y_zero), $"{nu},0");
            }

            for (int n = -15; n <= 15; n++) {
                foreach (ddouble alpha in new[] {
                    -double.ScaleB(1, -12), -double.ScaleB(1, -24), -double.ScaleB(1, -48), -double.ScaleB(1, -96),
                    double.ScaleB(1, -12), double.ScaleB(1, -24), double.ScaleB(1, -48), double.ScaleB(1, -96), }) {

                    ddouble nu = n + alpha;

                    for (ddouble x = ddouble.Ldexp(1, -950); x > 0; x /= 2) {
                        ddouble y = DDoubleOptimizedBessel.RealBessel.BesselI(nu, x);

                        Assert.IsFalse(ddouble.IsNaN(y), $"{nu},{x}");
                    }

                    ddouble y_zero = DDoubleOptimizedBessel.RealBessel.BesselI(nu, 0);

                    Assert.IsFalse(ddouble.IsNaN(y_zero), $"{nu},0");
                }
            }
        }

        [TestMethod]
        public void BesselKLimitZeroXTest() {
            ddouble eps = ddouble.Ldexp(1, -4);

            for (ddouble nu = -16; nu <= 16; nu += 1d / 32) {
                Console.WriteLine($"nu={nu}");

                ddouble y0 = DDoubleOptimizedBessel.RealBessel.BesselK(nu, eps * 4);
                ddouble y1 = DDoubleOptimizedBessel.RealBessel.BesselK(nu, eps * 2);

                if (ddouble.IsInfinity(y0) || ddouble.IsInfinity(y1)) {
                    continue;
                }

                for (ddouble x = eps; x > 0; x /= 2) {
                    ddouble y2 = DDoubleOptimizedBessel.RealBessel.BesselK(nu, x);

                    Console.WriteLine($"x = {x}\ny = {y2}");

                    Assert.IsFalse(ddouble.IsNaN(y2), $"{nu},{x}");

                    if (y1 > 1e-300) {
                        if (ddouble.Abs(y0) >= ddouble.Abs(y1)) {
                            Assert.IsTrue(ddouble.Abs(y1) >= ddouble.Abs(y2), $"{nu},{x}");
                        }
                        if (ddouble.Abs(y0) <= ddouble.Abs(y1)) {
                            Assert.IsTrue(ddouble.Abs(y1) <= ddouble.Abs(y2), $"{nu},{x}");
                        }
                    }

                    if (y2 == 0) {
                        break;
                    }

                    Assert.AreEqual(ddouble.Sign(y1), ddouble.Sign(y2), $"{nu},{x}");

                    if (ddouble.IsInfinity(y2)) {
                        break;
                    }

                    (y0, y1) = (y1, y2);
                }
            }

            for (ddouble nu = -16; nu <= 16; nu += 1d / 256) {
                for (ddouble x = ddouble.Ldexp(1, -950); x > 0; x /= 2) {
                    ddouble y = DDoubleOptimizedBessel.RealBessel.BesselK(nu, x);

                    Assert.IsFalse(ddouble.IsNaN(y), $"{nu},{x}");
                    Assert.IsTrue(ddouble.IsPositive(y), $"{nu},{x}");
                }

                ddouble y_zero = DDoubleOptimizedBessel.RealBessel.BesselK(nu, 0);

                Assert.IsFalse(ddouble.IsNaN(y_zero), $"{nu},0");
            }

            for (int n = -15; n <= 15; n++) {
                foreach (ddouble alpha in new[] {
                    -double.ScaleB(1, -12), -double.ScaleB(1, -24), -double.ScaleB(1, -48), -double.ScaleB(1, -96),
                    double.ScaleB(1, -12), double.ScaleB(1, -24), double.ScaleB(1, -48), double.ScaleB(1, -96), }) {

                    ddouble nu = n + alpha;

                    for (ddouble x = ddouble.Ldexp(1, -950); x > 0; x /= 2) {
                        ddouble y = DDoubleOptimizedBessel.RealBessel.BesselK(nu, x);

                        Assert.IsFalse(ddouble.IsNaN(y), $"{nu},{x}");
                        Assert.IsTrue(ddouble.IsPositive(y), $"{nu},{x}");
                    }

                    ddouble y_zero = DDoubleOptimizedBessel.RealBessel.BesselK(nu, 0);

                    Assert.IsFalse(ddouble.IsNaN(y_zero), $"{nu},0");
                }
            }
        }
    }
}
