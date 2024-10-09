using DoubleDouble;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DDoubleOptimizedBesselTests {
    [TestClass()]
    public class RealBesselAbnormalValueTests {
        [TestMethod]
        public void BesselJAbnormalTest() {
            for (ddouble nu = -16; nu <= 16; nu += 0.125) {
                Assert.IsTrue(ddouble.IsZero(DDoubleOptimizedBessel.RealBessel.BesselJ(nu, ddouble.PositiveInfinity)));
                Assert.IsTrue(ddouble.IsFinite(DDoubleOptimizedBessel.RealBessel.BesselJ(nu, ddouble.MaxValue)));
                Assert.IsTrue(ddouble.IsFinite(DDoubleOptimizedBessel.RealBessel.BesselJ(nu, 1e+300)));
                Assert.IsTrue(ddouble.IsNaN(DDoubleOptimizedBessel.RealBessel.BesselJ(nu, ddouble.NaN)));
            }
        }

        [TestMethod]
        public void BesselYAbnormalTest() {
            for (ddouble nu = -16; nu <= 16; nu += 0.125) {
                Assert.IsTrue(ddouble.IsZero(DDoubleOptimizedBessel.RealBessel.BesselY(nu, ddouble.PositiveInfinity)));
                Assert.IsTrue(ddouble.IsFinite(DDoubleOptimizedBessel.RealBessel.BesselY(nu, ddouble.MaxValue)));
                Assert.IsTrue(ddouble.IsFinite(DDoubleOptimizedBessel.RealBessel.BesselY(nu, 1e+300)));
                Assert.IsTrue(ddouble.IsNaN(DDoubleOptimizedBessel.RealBessel.BesselY(nu, ddouble.NaN)));
            }
        }

        [TestMethod]
        public void BesselIAbnormalTest() {
            for (ddouble nu = -16; nu <= 16; nu += 0.125) {
                Assert.IsTrue(ddouble.IsPositiveInfinity(DDoubleOptimizedBessel.RealBessel.BesselI(nu, ddouble.PositiveInfinity)));
                Assert.IsTrue(ddouble.IsPositiveInfinity(DDoubleOptimizedBessel.RealBessel.BesselI(nu, ddouble.MaxValue)));
                Assert.IsTrue(ddouble.IsPositiveInfinity(DDoubleOptimizedBessel.RealBessel.BesselI(nu, 1e+300)));
                Assert.IsTrue(ddouble.IsNaN(DDoubleOptimizedBessel.RealBessel.BesselI(nu, ddouble.NaN)));

                Assert.IsTrue(ddouble.IsPlusZero(DDoubleOptimizedBessel.RealBessel.BesselI(nu, ddouble.PositiveInfinity, scale: true)));
                Assert.IsTrue(ddouble.IsFinite(DDoubleOptimizedBessel.RealBessel.BesselI(nu, ddouble.MaxValue, scale: true)));
                Assert.IsTrue(ddouble.IsFinite(DDoubleOptimizedBessel.RealBessel.BesselI(nu, 1e+300, scale: true)));
                Assert.IsTrue(ddouble.IsNaN(DDoubleOptimizedBessel.RealBessel.BesselI(nu, ddouble.NaN, scale: true)));
            }
        }

        [TestMethod]
        public void BesselKAbnormalTest() {
            for (ddouble nu = 0; nu <= 16; nu += 0.125) {
                Assert.IsTrue(ddouble.IsPlusZero(DDoubleOptimizedBessel.RealBessel.BesselK(nu, ddouble.PositiveInfinity)));
                Assert.IsTrue(ddouble.IsFinite(DDoubleOptimizedBessel.RealBessel.BesselK(nu, ddouble.MaxValue)));
                Assert.IsTrue(ddouble.IsFinite(DDoubleOptimizedBessel.RealBessel.BesselK(nu, 1e+300)));
                Assert.IsTrue(ddouble.IsNaN(DDoubleOptimizedBessel.RealBessel.BesselK(nu, ddouble.NaN)));

                Assert.IsTrue(ddouble.IsPlusZero(DDoubleOptimizedBessel.RealBessel.BesselK(nu, ddouble.PositiveInfinity, scale: true)));
                Assert.IsTrue(ddouble.IsFinite(DDoubleOptimizedBessel.RealBessel.BesselK(nu, ddouble.MaxValue, scale: true)));
                Assert.IsTrue(ddouble.IsFinite(DDoubleOptimizedBessel.RealBessel.BesselK(nu, 1e+300, scale: true)));
                Assert.IsTrue(ddouble.IsNaN(DDoubleOptimizedBessel.RealBessel.BesselK(nu, ddouble.NaN, scale: true)));
            }
        }

        [TestMethod]
        public void BesselJZeroTest() {
            ddouble[] x0_expecteds = {
                ddouble.Zero,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.Zero,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.Zero,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.Zero,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.Zero,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.Zero,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.Zero,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.Zero,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.Zero,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.Zero,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.Zero,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.Zero,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.Zero,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.Zero,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.Zero,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.Zero,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                1,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
            };

            foreach ((ddouble x, ddouble[] expecteds) in new (ddouble, ddouble[])[] {(0, x0_expecteds)}) {

                for ((ddouble nu, int i) = (-16, 0); i < expecteds.Length; nu += 0.25d, i++) {
                    ddouble expected = expecteds[i];

                    ddouble actual = DDoubleOptimizedBessel.RealBessel.BesselJ(nu, x);
                    Assert.AreEqual(expected, actual, $"{x},{nu}");
                }
            }
        }

        [TestMethod]
        public void BesselYZeroTest() {
            ddouble[] x0_expecteds = {
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.Zero,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.Zero,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.Zero,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.Zero,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.Zero,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.Zero,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.Zero,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.Zero,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.Zero,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.Zero,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.Zero,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.Zero,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.Zero,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.Zero,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.Zero,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.Zero,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
            };

            foreach ((ddouble x, ddouble[] expecteds) in new (ddouble, ddouble[])[] {(0, x0_expecteds)}) {

                for ((ddouble nu, int i) = (-16, 0); i < expecteds.Length; nu += 0.25d, i++) {
                    ddouble expected = expecteds[i];

                    ddouble actual = DDoubleOptimizedBessel.RealBessel.BesselY(nu, x);
                    Assert.AreEqual(expected, actual, $"{x},{nu}");
                }
            }
        }

        [TestMethod]
        public void BesselIZeroTest() {
            ddouble[] x0_expecteds = {
                ddouble.Zero,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.Zero,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.Zero,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.Zero,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.Zero,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.Zero,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.Zero,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.Zero,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.Zero,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.Zero,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.Zero,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.Zero,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.Zero,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.Zero,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.Zero,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.NegativeInfinity,
                ddouble.Zero,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                1,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
                ddouble.Zero,
            };

            foreach ((ddouble x, ddouble[] expecteds) in new (ddouble, ddouble[])[] {(0, x0_expecteds)}) {

                for ((ddouble nu, int i) = (-16, 0); i < expecteds.Length; nu += 0.25d, i++) {
                    ddouble expected = expecteds[i];

                    ddouble actual = DDoubleOptimizedBessel.RealBessel.BesselI(nu, x);
                    Assert.AreEqual(expected, actual, $"{x},{nu}");
                }
            }
        }

        [TestMethod]
        public void BesselKZeroTest() {
            ddouble[] x0_expecteds = {
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
                ddouble.PositiveInfinity,
            };

            foreach ((ddouble x, ddouble[] expecteds) in new (ddouble, ddouble[])[] {(0, x0_expecteds)}) {

                for ((ddouble nu, int i) = (0, 0); i < expecteds.Length; nu += 0.25d, i++) {
                    ddouble expected = expecteds[i];

                    ddouble actual = DDoubleOptimizedBessel.RealBessel.BesselK(nu, x);
                    Assert.AreEqual(expected, actual, $"{x},{nu}");
                }
            }
        }
    }
}
