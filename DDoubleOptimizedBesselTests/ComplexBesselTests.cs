﻿using ComplexBessel;
using DoubleDouble;
using DoubleDoubleComplex;

namespace DDoubleOptimizedBesselTests {
    [TestClass()]
    public class ComplexBesselTests {
        readonly Complex[] zs_mini = [
            (64, -8), (64, 8), (-64, -8), (-64, 8), (-8, 64), (8, 64), (-8, -64), (8, -64)
        ];

        [TestMethod()]
        public void BesselJTest() {
            for (double nu = -16; nu <= 16; nu += 0.25) {
                Console.WriteLine(nu);

                for (double r = 0; r <= 42; r += 0.5) {
                    for (double i = 0; i <= 42; i += 0.5) {
                        if (r == 0 && i == 0) {
                            continue;
                        }

                        Complex expected = BesselN4.BesselJ(nu, (r, i)).ToString();
                        Complex actual = DDoubleOptimizedBessel.ComplexBessel.BesselJ(nu, (r, i));

                        ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                        Console.WriteLine($"{nu}, {(r, i)}, {err:e4}");
                        Console.WriteLine(expected);
                        Console.WriteLine(actual);

                        Assert.IsTrue(err < 2e-27, $"\n{nu}, {(r, i)}\n{expected}\n{actual}\n{err:e4}");
                    }
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselYTest() {
            for (double nu = -16; nu <= 16; nu += 0.25) {
                Console.WriteLine(nu);

                for (double r = 0; r <= 42; r += 0.5) {
                    for (double i = 0; i <= 42; i += 0.5) {
                        if (r == 0 && i == 0) {
                            continue;
                        }

                        Complex expected = BesselN4.BesselY(nu, (r, i)).ToString();
                        Complex actual = DDoubleOptimizedBessel.ComplexBessel.BesselY(nu, (r, i));

                        ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                        Console.WriteLine($"{nu}, {(r, i)}, {err:e4}");
                        Console.WriteLine(expected);
                        Console.WriteLine(actual);

                        Assert.IsTrue(err < 2e-27, $"\n{nu}, {(r, i)}\n{expected}\n{actual}\n{err:e4}");
                    }
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselITest() {
            for (double nu = -16; nu <= 16; nu += 0.25) {
                Console.WriteLine(nu);

                for (double r = 0; r <= 42; r += 0.5) {
                    for (double i = 0; i <= 42; i += 0.5) {
                        if (r == 0 && i == 0) {
                            continue;
                        }

                        Complex expected = BesselN4.BesselI(nu, (r, i)).ToString();
                        Complex actual = DDoubleOptimizedBessel.ComplexBessel.BesselI(nu, (r, i));

                        ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                        Console.WriteLine($"{nu}, {(r, i)}, {err:e4}");
                        Console.WriteLine(expected);
                        Console.WriteLine(actual);

                        Assert.IsTrue(err < 2e-27, $"\n{nu}, {(r, i)}\n{expected}\n{actual}\n{err:e4}");
                    }
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselKTest() {
            for (double nu = 0; nu <= 16; nu += 0.25) {
                Console.WriteLine(nu);

                for (double r = 0; r <= 42; r += 0.5) {
                    for (double i = 0; i <= 42; i += 0.5) {
                        if (r == 0 && i == 0) {
                            continue;
                        }

                        Complex expected = BesselN4.BesselK(nu, (r, i)).ToString();
                        Complex actual = DDoubleOptimizedBessel.ComplexBessel.BesselK(nu, (r, i));

                        ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                        Console.WriteLine($"{nu}, {(r, i)}, {err:e4}");
                        Console.WriteLine(expected);
                        Console.WriteLine(actual);

                        Assert.IsTrue(err < 2e-27, $"\n{nu}, {(r, i)}\n{expected}\n{actual}\n{err:e4}");
                    }
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselYNearIntTest() {
            for (int n = -16; n <= 16; n++) {
                for (double alpha = 0; alpha <= 0.25; alpha = alpha < double.ScaleB(1, -8) ? double.ScaleB(1, -8) : alpha * 2) {
                    Console.WriteLine(n + alpha);

                    for (double r = 1d / 4; r <= 8; r += 1d / 4) {
                        for (double i = 1d / 4; i <= 8; i += 1d / 4) {
                            Complex expected = BesselN4.BesselY(n + alpha, (r, i)).ToString();
                            Complex actual = DDoubleOptimizedBessel.ComplexBessel.BesselY(n + alpha, (r, i));

                            ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                            Console.WriteLine($"{n + alpha}, {(r, i)}, {err:e4}");
                            Console.WriteLine(expected);
                            Console.WriteLine(actual);

                            Assert.IsTrue(err < 4e-27, $"\n{n + alpha}, {(r, i)}\n{expected}\n{actual}\n{err:e4}");
                        }
                    }

                    if (alpha == 0d) {
                        continue;
                    }

                    Console.WriteLine(n - alpha);

                    for (double r = 1d / 4; r <= 8; r += 1d / 4) {
                        for (double i = 1d / 4; i <= 8; i += 1d / 4) {
                            Complex expected = BesselN4.BesselY(n - alpha, (r, i)).ToString();
                            Complex actual = DDoubleOptimizedBessel.ComplexBessel.BesselY(n - alpha, (r, i));

                            ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                            Console.WriteLine($"{n - alpha}, {(r, i)}, {err:e4}");
                            Console.WriteLine(expected);
                            Console.WriteLine(actual);

                            Assert.IsTrue(err < 4e-27, $"\n{n - alpha}, {(r, i)}\n{expected}\n{actual}\n{err:e4}");
                        }
                    }

                    Console.WriteLine(string.Empty);
                }
            }
        }

        [TestMethod()]
        public void BesselKNearIntTest() {
            for (int n = 0; n <= 16; n++) {
                for (double alpha = 0; alpha <= 0.25; alpha = alpha < double.ScaleB(1, -8) ? double.ScaleB(1, -8) : alpha * 2) {
                    Console.WriteLine(n + alpha);

                    for (double r = 1d / 8; r <= 4; r += 1d / 8) {
                        for (double i = 1d / 8; i <= 4; i += 1d / 8) {
                            Complex expected = BesselN4.BesselK(n + alpha, (r, i)).ToString();
                            Complex actual = DDoubleOptimizedBessel.ComplexBessel.BesselK(n + alpha, (r, i));

                            ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                            Console.WriteLine($"{n + alpha}, {(r, i)}, {err:e4}");
                            Console.WriteLine(expected);
                            Console.WriteLine(actual);

                            Assert.IsTrue(err < 4e-27, $"\n{n + alpha}, {(r, i)}\n{expected}\n{actual}\n{err:e4}");
                        }
                    }

                    if (n <= 0 || alpha == 0d) {
                        continue;
                    }

                    Console.WriteLine(n - alpha);

                    for (double r = 1d / 8; r <= 4; r += 1d / 8) {
                        for (double i = 1d / 8; i <= 4; i += 1d / 8) {
                            Complex expected = BesselN4.BesselK(n - alpha, (r, i)).ToString();
                            Complex actual = DDoubleOptimizedBessel.ComplexBessel.BesselK(n - alpha, (r, i));

                            ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                            Console.WriteLine($"{n - alpha}, {(r, i)}, {err:e4}");
                            Console.WriteLine(expected);
                            Console.WriteLine(actual);

                            Assert.IsTrue(err < 4e-27, $"\n{n - alpha}, {(r, i)}\n{expected}\n{actual}\n{err:e4}");
                        }
                    }

                    Console.WriteLine(string.Empty);
                }
            }
        }

        [TestMethod()]
        public void BesselJMinusRTest() {
            for (double nu = -16; nu <= 16; nu += 0.25) {
                Console.WriteLine(nu);

                foreach (double r in new double[] { 0, -1 / 8d, -1 / 4d, -1 / 2d, -1, -2, -4, -8, -16, -32, -64 }) {
                    foreach (double i in new double[] { 1 / 8d, 1 / 4d, 1 / 2d, 1, 2, 4, 8, 16, 32, 64 }) {

                        Complex expected = BesselN4.BesselJ(nu, (r, i)).ToString();
                        Complex actual = DDoubleOptimizedBessel.ComplexBessel.BesselJ(nu, (r, i));

                        ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                        Console.WriteLine($"{nu}, {(r, i)}, {err:e4}");
                        Console.WriteLine(expected);
                        Console.WriteLine(actual);

                        Assert.IsTrue(err < 4e-27, $"\n{nu}, {(r, i)}\n{expected}\n{actual}\n{err:e4}");
                    }
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselYMinusRTest() {
            for (double nu = -16; nu <= 16; nu += 0.25) {
                Console.WriteLine(nu);

                foreach (double r in new double[] { 0, -1 / 8d, -1 / 4d, -1 / 2d, -1, -2, -4, -8, -16, -32, -64 }) {
                    foreach (double i in new double[] { 1 / 8d, 1 / 4d, 1 / 2d, 1, 2, 4, 8, 16, 32, 64 }) {

                        Complex expected = BesselN4.BesselY(nu, (r, i)).ToString();
                        Complex actual = DDoubleOptimizedBessel.ComplexBessel.BesselY(nu, (r, i));

                        ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                        Console.WriteLine($"{nu}, {(r, i)}, {err:e4}");
                        Console.WriteLine(expected);
                        Console.WriteLine(actual);

                        Assert.IsTrue(err < 2e-27, $"\n{nu}, {(r, i)}\n{expected}\n{actual}\n{err:e4}");
                    }
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselIMinusRTest() {
            for (double nu = -16; nu <= 16; nu += 0.25) {
                Console.WriteLine(nu);

                foreach (double r in new double[] { 0, -1 / 8d, -1 / 4d, -1 / 2d, -1, -2, -4, -8, -16, -32, -64 }) {
                    foreach (double i in new double[] { 1 / 8d, 1 / 4d, 1 / 2d, 1, 2, 4, 8, 16, 32, 64 }) {

                        Complex expected = BesselN4.BesselI(nu, (r, i)).ToString();
                        Complex actual = DDoubleOptimizedBessel.ComplexBessel.BesselI(nu, (r, i));

                        ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                        Console.WriteLine($"{nu}, {(r, i)}, {err:e4}");
                        Console.WriteLine(expected);
                        Console.WriteLine(actual);

                        Assert.IsTrue(err < 4e-27, $"\n{nu}, {(r, i)}\n{expected}\n{actual}\n{err:e4}");
                    }
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselKMinusRTest() {
            for (double nu = 0; nu <= 16; nu += 0.25) {
                Console.WriteLine(nu);

                foreach (double r in new double[] { 0, -1 / 8d, -1 / 4d, -1 / 2d, -1, -2, -4, -8, -16, -32, -64 }) {
                    foreach (double i in new double[] { 1 / 8d, 1 / 4d, 1 / 2d, 1, 2, 4, 8, 16, 32, 64 }) {

                        Complex expected = BesselN4.BesselK(nu, (r, i)).ToString();
                        Complex actual = DDoubleOptimizedBessel.ComplexBessel.BesselK(nu, (r, i));

                        ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                        Console.WriteLine($"{nu}, {(r, i)}, {err:e4}");
                        Console.WriteLine(expected);
                        Console.WriteLine(actual);

                        Assert.IsTrue(err < 2e-27, $"\n{nu}, {(r, i)}\n{expected}\n{actual}\n{err:e4}");
                    }
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselJRecurrenceTest() {
            double[] vs = [0, 1 / 8d, 1 / 4d, 1 / 2d, 1, 2, 4, 8, 16, 32, 64, 128, 256];

            double[] nus = [16.25, 16.5, 16.75, 17, 17.25, 17.5, 17.75, 18,
                           18.25, 18.5, 18.75, 19, 19.5, 20, 20.25, 20.5, 20.75,
                           21, 63.75, 64, 127.75, 128, 255.5, 255.75, 256];

            foreach (double nu_abs in nus) {

                foreach (double nu in new double[] { nu_abs, -nu_abs }) {
                    Console.WriteLine(nu);

                    foreach (double r in vs) {
                        foreach (double i in vs) {
                            if (r == 0 && i == 0) {
                                continue;
                            }

                            Complex expected = BesselN4.BesselJ(nu, (r, i)).ToString();
                            Complex actual = DDoubleOptimizedBessel.ComplexBessel.BesselJ(nu, (r, i));

                            ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                            Console.WriteLine($"{nu}, {(r, i)}, {err:e4}");
                            Console.WriteLine(expected);
                            Console.WriteLine(actual);

                            if (!Complex.IsFinite(expected)) {
                                continue; //ignore
                            }

                            if (Complex.IsZero(expected)) {
                                Assert.IsTrue(Complex.IsZero(actual));
                                continue;
                            }

                            if (expected.Magnitude < 1e-292) {
                                Assert.IsTrue(actual.Magnitude < 1e-291);
                                continue;
                            }

                            Assert.IsTrue(err < 4e-27, $"\n{nu}, {(r, i)}\n{expected}\n{actual}\n{err:e4}");
                        }
                    }

                    Console.WriteLine(string.Empty);
                }
            }
        }

        [TestMethod()]
        public void BesselYRecurrenceTest() {
            double[] vs = [0, 1 / 8d, 1 / 4d, 1 / 2d, 1, 2, 4, 8, 16, 32, 64, 128, 256];

            double[] nus = [16.25, 16.5, 16.75, 17, 17.25, 17.5, 17.75, 18,
                           18.25, 18.5, 18.75, 19, 19.5, 20, 20.25, 20.5, 20.75,
                           21, 63.75, 64, 127.75, 128, 255.5, 255.75, 256];

            foreach (double nu_abs in nus) {

                foreach (double nu in new double[] { nu_abs, -nu_abs }) {
                    Console.WriteLine(nu);

                    foreach (double r in vs) {
                        foreach (double i in vs) {
                            if (r == 0 && i == 0) {
                                continue;
                            }

                            Complex expected = BesselN4.BesselY(nu, (r, i)).ToString();
                            Complex actual = DDoubleOptimizedBessel.ComplexBessel.BesselY(nu, (r, i));

                            ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                            Console.WriteLine($"{nu}, {(r, i)}, {err:e4}");
                            Console.WriteLine(expected);
                            Console.WriteLine(actual);

                            if (!Complex.IsFinite(expected)) {
                                continue; //ignore
                            }

                            if (Complex.IsZero(expected)) {
                                Assert.IsTrue(Complex.IsZero(actual));
                                continue;
                            }

                            if (expected.Magnitude < 1e-292) {
                                Assert.IsTrue(actual.Magnitude < 1e-291);
                                continue;
                            }

                            Assert.IsTrue(err < 4e-27, $"\n{nu}, {(r, i)}\n{expected}\n{actual}\n{err:e4}");
                        }
                    }

                    Console.WriteLine(string.Empty);
                }
            }
        }

        [TestMethod()]
        public void BesselIRecurrenceTest() {
            double[] vs = [0, 1 / 8d, 1 / 4d, 1 / 2d, 1, 2, 4, 8, 16, 32, 64, 128, 256];

            double[] nus = [16.25, 16.5, 16.75, 17, 17.25, 17.5, 17.75, 18,
                           18.25, 18.5, 18.75, 19, 19.5, 20, 20.25, 20.5, 20.75,
                           21, 63.75, 64, 127.75, 128, 255.5, 255.75, 256];

            foreach (double nu_abs in nus) {

                foreach (double nu in new double[] { nu_abs, -nu_abs }) {
                    Console.WriteLine(nu);

                    foreach (double r in vs) {
                        foreach (double i in vs) {
                            if (r == 0 && i == 0) {
                                continue;
                            }

                            Complex expected = BesselN4.BesselI(nu, (r, i)).ToString();
                            Complex actual = DDoubleOptimizedBessel.ComplexBessel.BesselI(nu, (r, i));

                            ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                            Console.WriteLine($"{nu}, {(r, i)}, {err:e4}");
                            Console.WriteLine(expected);
                            Console.WriteLine(actual);

                            if (!Complex.IsFinite(expected)) {
                                continue; //ignore
                            }

                            if (Complex.IsZero(expected)) {
                                Assert.IsTrue(Complex.IsZero(actual));
                                continue;
                            }

                            if (expected.Magnitude < 1e-292) {
                                Assert.IsTrue(actual.Magnitude < 1e-291);
                                continue;
                            }

                            Assert.IsTrue(err < 4e-27, $"\n{nu}, {(r, i)}\n{expected}\n{actual}\n{err:e4}");
                        }
                    }

                    Console.WriteLine(string.Empty);
                }
            }
        }

        [TestMethod()]
        public void BesselKRecurrenceTest() {
            double[] vs = [0, 1 / 8d, 1 / 4d, 1 / 2d, 1, 2, 4, 8, 16, 32, 64, 128, 256];

            double[] nus = [16.25, 16.5, 16.75, 17, 17.25, 17.5, 17.75, 18,
                           18.25, 18.5, 18.75, 19, 19.5, 20, 20.25, 20.5, 20.75,
                           21, 63.75, 64, 127.75, 128, 255.5, 255.75, 256];

            foreach (double nu_abs in nus) {

                foreach (double nu in new double[] { nu_abs, -nu_abs }) {
                    Console.WriteLine(nu);

                    foreach (double r in vs) {
                        foreach (double i in vs) {
                            if (r == 0 && i == 0) {
                                continue;
                            }

                            Complex expected = BesselN4.BesselK(nu, (r, i)).ToString();

                            Complex actual = DDoubleOptimizedBessel.ComplexBessel.BesselK(nu, (r, i));

                            ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                            Console.WriteLine($"{nu}, {(r, i)}, {err:e4}");
                            Console.WriteLine(expected);
                            Console.WriteLine(actual);

                            if (!Complex.IsFinite(expected)) {
                                continue; //ignore
                            }

                            if (Complex.IsZero(expected)) {
                                Assert.IsTrue(Complex.IsZero(actual));
                                continue;
                            }

                            if (expected.Magnitude < 1e-292) {
                                Assert.IsTrue(actual.Magnitude < 1e-291);
                                continue;
                            }

                            Assert.IsTrue(err < 2e-27, $"\n{nu}, {(r, i)}\n{expected}\n{actual}\n{err:e4}");
                        }
                    }

                    Console.WriteLine(string.Empty);
                }
            }
        }

        [TestMethod()]
        public void HankelH1Test() {
            for (double nu = -16; nu <= 16; nu += 0.25) {
                Console.WriteLine(nu);

                double[] xs = [
                    0, 1 / 8, 1 / 4, 1 / 2, 1, 2, 4, 8, 16, 32, 64,
                    -1 / 8, -1 / 4, -1 / 2, -1, -2, -4, -8, -16, -32, -64
                ];

                foreach (double r in xs) {
                    foreach (double i in xs) {
                        if (r == 0 && i <= 0) {
                            continue;
                        }

                        Complex expected = BesselN4.HankelH1(nu, (r, i)).ToString();
                        Complex actual = DDoubleOptimizedBessel.ComplexBessel.HankelH1(nu, (r, i));

                        ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                        Console.WriteLine($"{nu}, {(r, i)}, {err:e4}");
                        Console.WriteLine(expected);
                        Console.WriteLine(actual);

                        if (double.Abs(nu) != 1.5 || (r, i) != (0, -1)) {
                            Assert.IsTrue(err < 2e-27, $"\n{nu}, {(r, i)}\n{expected}\n{actual}\n{err}");
                        }
                        else {
                            Assert.IsTrue(actual.Magnitude < 1e-30, $"\n{nu}, {(r, i)}\n{expected}\n{actual}\n{err}");
                        }
                    }
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void HankelH2Test() {
            for (double nu = -16; nu <= 16; nu += 0.25) {
                Console.WriteLine(nu);

                double[] xs = [
                    0, 1 / 8, 1 / 4, 1 / 2, 1, 2, 4, 8, 16, 32, 64,
                    -1 / 8, -1 / 4, -1 / 2, -1, -2, -4, -8, -16, -32, -64
                ];

                foreach (double r in xs) {
                    foreach (double i in xs) {
                        if (r == 0 && i <= 0) {
                            continue;
                        }

                        Complex expected = BesselN4.HankelH2(nu, (r, i)).ToString();
                        Complex actual = DDoubleOptimizedBessel.ComplexBessel.HankelH2(nu, (r, i));

                        ddouble err = (expected - actual).Magnitude / expected.Magnitude;

                        Console.WriteLine($"{nu}, {(r, i)}, {err:e4}");
                        Console.WriteLine(expected);
                        Console.WriteLine(actual);

                        if (double.Abs(nu) != 1.5 || (r, i) != (0, 1)) {
                            Assert.IsTrue(err < 2e-27, $"\n{nu}, {(r, i)}\n{expected}\n{actual}\n{err}");
                        }
                        else {
                            Assert.IsTrue(actual.Magnitude < 1e-30, $"\n{nu}, {(r, i)}\n{expected}\n{actual}\n{err}");
                        }
                    }
                }

                Console.WriteLine(string.Empty);
            }
        }

        [TestMethod()]
        public void BesselJNu1p25Test() {
            Complex[] expecteds = [
                "9.2145929015511005131094720273851178236-147.6159570273588415961403709192966375974i",
                "9.2145929015511005131094720273851178236+147.6159570273588415961403709192966375974i",
                "-110.8959453519476388434291742350541699970-97.8645430988272226734726977555188560612i",
                "-110.8959453519476388434291742350541699970+97.8645430988272226734726977555188560612i",
                "-2.722560899830383623022272861089442378455e26-1.409865452209600514674146268526138161281e26i",
                "2.922066696281459872639446880884572748286e26+9.28215852645366138489133371390795409174e25i",
                "-2.722560899830383623022272861089442378455e26+1.409865452209600514674146268526138161281e26i",
                "2.922066696281459872639446880884572748286e26-9.28215852645366138489133371390795409174e25i"
            ];

            foreach ((Complex z, Complex expected) in zs_mini.Zip(expecteds)) {
                Complex actual = DDoubleOptimizedBessel.ComplexBessel.BesselJ(1.25, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 2e-30);
            }
        }

        [TestMethod()]
        public void BesselYNu1p25Test() {
            Complex[] expecteds = [
                "-147.6159903358803008197446577402684011815-9.2145949710669354056471653366179386559i",
                "-147.6159903358803008197446577402684011815+9.2145949710669354056471653366179386559i",
                "-97.8645180827771469333720266685710849958+110.8959674412603533326241154327892910921i",
                "-97.8645180827771469333720266685710849958-110.8959674412603533326241154327892910921i",
                "1.409865452209600514674146268526138161281e26-2.722560899830383623022272861089442378455e26i",
                "-9.28215852645366138489133371390795409174e25+2.922066696281459872639446880884572748286e26i",
                "1.409865452209600514674146268526138161281e26+2.722560899830383623022272861089442378455e26i",
                "-9.28215852645366138489133371390795409174e25-2.922066696281459872639446880884572748286e26i"
            ];

            foreach ((Complex z, Complex expected) in zs_mini.Zip(expecteds)) {
                Complex actual = DDoubleOptimizedBessel.ComplexBessel.BesselY(1.25, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 2e-30);
            }
        }

        [TestMethod()]
        public void BesselINu1p25Test() {
            Complex[] expecteds = [
                "-2.60666884921141062870887160247972666843e25-3.054850441793332396273499925212507250408e26i",
                "-2.60666884921141062870887160247972666843e25+3.054850441793332396273499925212507250408e26i",
                "-1.975786140944273667434472492382389247086e26-2.344424784861298183212213072859469902627e26i",
                "-1.975786140944273667434472492382389247086e26+2.344424784861298183212213072859469902627e26i",
                "-139.9056334090550484779115612126481548920-47.9770073249204233223395686709723050197i",
                "132.8530893302299147044482156602733629886+65.0033548892541275819643373472106188313i",
                "-139.9056334090550484779115612126481548920+47.9770073249204233223395686709723050197i",
                "132.8530893302299147044482156602733629886-65.0033548892541275819643373472106188313i"
            ];

            foreach ((Complex z, Complex expected) in zs_mini.Zip(expecteds)) {
                Complex actual = DDoubleOptimizedBessel.ComplexBessel.BesselI(1.25, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 2e-30);
            }
        }

        [TestMethod()]
        public void BesselKNu1p25Test() {
            Complex[] expecteds = [
                "-5.25616659269784449308681704247624052439e-30+2.472840504427643132508872150883116217499e-29i",
                "-5.25616659269784449308681704247624052439e-30-2.472840504427643132508872150883116217499e-29i",
                "-9.597095705753467335604038328256260292002e26-8.18909170702392813709399161901292238171e25i",
                "-9.597095705753467335604038328256260292002e26+8.18909170702392813709399161901292238171e25i",
                "-204.2140852044501196181755107589505155153-417.3703365407477352155598375752493218656i",
                "-0.00001701900642617936969720732433114038442166-0.00004958223421184402515747294308650095460414i",
                "-204.2140852044501196181755107589505155153+417.3703365407477352155598375752493218656i",
                "-0.00001701900642617936969720732433114038442166+0.00004958223421184402515747294308650095460414i"
            ];

            foreach ((Complex z, Complex expected) in zs_mini.Zip(expecteds)) {
                Complex actual = DDoubleOptimizedBessel.ComplexBessel.BesselK(1.25, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 2e-30);
            }
        }

        [TestMethod()]
        public void BesselJNuM1p75Test() {
            Complex[] expecteds = [
                "-147.2999139438063945323394630815032224390-10.9056339300812418071382865800802817172i",
                "-147.2999139438063945323394630815032224390+10.9056339300812418071382865800802817172i",
                "-96.4453203127618438412463729175939687224+111.8682157229589337900597395166247943914i",
                "-96.4453203127618438412463729175939687224-111.8682157229589337900597395166247943914i",
                "1.397508963018419728866125463198502317718e26-2.689047103063237008011487557828466479775e26i",
                "-9.13255376986751184684287116150543446463e25+2.889631506025360495813279339858706432795e26i",
                "1.397508963018419728866125463198502317718e26+2.689047103063237008011487557828466479775e26i",
                "-9.13255376986751184684287116150543446463e25-2.889631506025360495813279339858706432795e26i"
            ];

            foreach ((Complex z, Complex expected) in zs_mini.Zip(expecteds)) {
                Complex actual = DDoubleOptimizedBessel.ComplexBessel.BesselJ(-1.75, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 2e-30);
            }
        }

        [TestMethod()]
        public void BesselYNuM1p75Test() {
            Complex[] expecteds = [
                "-10.9056322416899217958437967766265351684+147.2998805626433483749635926306468094384i",
                "-10.9056322416899217958437967766265351684-147.2998805626433483749635926306468094384i",
                "111.8681933127851316348824917693777891420+96.4453451106815493494087232567175383028i",
                "111.8681933127851316348824917693777891420-96.4453451106815493494087232567175383028i",
                "2.689047103063237008011487557828466479775e26+1.397508963018419728866125463198502317718e26i",
                "-2.889631506025360495813279339858706432795e26-9.13255376986751184684287116150543446463e25i",
                "2.689047103063237008011487557828466479775e26-1.397508963018419728866125463198502317718e26i",
                "-2.889631506025360495813279339858706432795e26+9.13255376986751184684287116150543446463e25i"
            ];

            foreach ((Complex z, Complex expected) in zs_mini.Zip(expecteds)) {
                Complex actual = DDoubleOptimizedBessel.ComplexBessel.BesselY(-1.75, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 2e-30);
            }
        }

        [TestMethod()]
        public void BesselINuM1p75Test() {
            Complex[] expecteds = [
                "-2.62076152242149906566470469518604540524e25-3.019159107207759513123694167620311528728e26i",
                "-2.62076152242149906566470469518604540524e25+3.019159107207759513123694167620311528728e26i",
                "1.949552053750027174411954077068820097967e26+2.320183702625431595125833192684946828683e26i",
                "1.949552053750027174411954077068820097967e26-2.320183702625431595125833192684946828683e26i",
                "140.2607810578372943130995226449327728200+46.2937446780354935621573347839373412787i",
                "131.9139702088759482017693824281792414315+66.4447286321608627589382194728601568628i",
                "140.2607810578372943130995226449327728200-46.2937446780354935621573347839373412787i",
                "131.9139702088759482017693824281792414315-66.4447286321608627589382194728601568628i"
            ];

            foreach ((Complex z, Complex expected) in zs_mini.Zip(expecteds)) {
                Complex actual = DDoubleOptimizedBessel.ComplexBessel.BesselI(-1.75, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 2e-30);
            }
        }

        [TestMethod()]
        public void BesselKNuM1p75Test() {
            Complex[] expecteds = [
                "-5.35222404565129229070775331686841841211e-30+2.500563373120195349488734686337018462950e-29i",
                "-5.35222404565129229070775331686841841211e-30-2.500563373120195349488734686337018462950e-29i",
                "-9.484968071222616258294302283203169493183e26-8.23336514565018365948586434406526860170e25i",
                "-9.484968071222616258294302283203169493183e26+8.23336514565018365948586434406526860170e25i",
                "-208.7422488243166520251860480998194073164-414.4199122853585266887857775703335245101i",
                "-0.00001761577059679953277111031956241000418243-0.00004945855291014150741401368246876143340060i",
                "-208.7422488243166520251860480998194073164+414.4199122853585266887857775703335245101i",
                "-0.00001761577059679953277111031956241000418243+0.00004945855291014150741401368246876143340060i"
            ];

            foreach ((Complex z, Complex expected) in zs_mini.Zip(expecteds)) {
                Complex actual = DDoubleOptimizedBessel.ComplexBessel.BesselK(-1.75, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 2e-30);
            }
        }

        [TestMethod()]
        public void BesselJNu1p375Test() {
            Complex[] expecteds = [
                "-19.3844491315920837390199087861121994352-146.5832933973515097634608469293879648659i",
                "-19.3844491315920837390199087861121994352+146.5832933973515097634608469293879648659i",
                "-128.0071970497257571256627952069408871441-74.0038936463615894478416609665561282765i",
                "-128.0071970497257571256627952069408871441+74.0038936463615894478416609665561282765i",
                "-2.388500783627850265566881132161440929619e26-1.909824495175187236085662122110732590600e26i",
                "2.678487439866468568876646704464127999900e26+1.475828794352374901455927115895960295955e26i",
                "-2.388500783627850265566881132161440929619e26+1.909824495175187236085662122110732590600e26i",
                "2.678487439866468568876646704464127999900e26-1.475828794352374901455927115895960295955e26i"
            ];

            foreach ((Complex z, Complex expected) in zs_mini.Zip(expecteds)) {
                Complex actual = DDoubleOptimizedBessel.ComplexBessel.BesselJ(1.375, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 2e-30);
            }
        }

        [TestMethod()]
        public void BesselYNu1p375Test() {
            Complex[] expecteds = [
                "-146.5833256942998126925823267036391419176+19.3844406821480655811156969072786874160i",
                "-146.5833256942998126925823267036391419176-19.3844406821480655811156969072786874160i",
                "-74.0038734805861684890514483972788984304+128.0072236547530183316730182054670692015i",
                "-74.0038734805861684890514483972788984304-128.0072236547530183316730182054670692015i",
                "1.909824495175187236085662122110732590600e26-2.388500783627850265566881132161440929619e26i",
                "-1.475828794352374901455927115895960295955e26+2.678487439866468568876646704464127999900e26i",
                "1.909824495175187236085662122110732590600e26+2.388500783627850265566881132161440929619e26i",
                "-1.475828794352374901455927115895960295955e26-2.678487439866468568876646704464127999900e26i"
            ];

            foreach ((Complex z, Complex expected) in zs_mini.Zip(expecteds)) {
                Complex actual = DDoubleOptimizedBessel.ComplexBessel.BesselY(1.375, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 2e-30);
            }
        }

        [TestMethod()]
        public void BesselINu1p375Test() {
            Complex[] expecteds = [
                "-2.60981095641589821397947573363886292337e25-3.047007460358397106095076131919830300429e26i",
                "-2.60981095641589821397947573363886292337e25+3.047007460358397106095076131919830300429e26i",
                "-2.715194686571793667532226846517769785438e26-1.407154366007622034466960613285609780914e26i",
                "-2.715194686571793667532226846517769785438e26+1.407154366007622034466960613285609780914e26i",
                "-111.1101312101308589510940681955478744057-97.5548948736905780898205433405748677084i",
                "132.6489770521213388555665512843130872017+65.3197340654040189559932767827654910729i",
                "-111.1101312101308589510940681955478744057+97.5548948736905780898205433405748677084i",
                "132.6489770521213388555665512843130872017-65.3197340654040189559932767827654910729i"
            ];

            foreach ((Complex z, Complex expected) in zs_mini.Zip(expecteds)) {
                Complex actual = DDoubleOptimizedBessel.ComplexBessel.BesselI(1.375, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 2e-30);
            }
        }

        [TestMethod()]
        public void BesselKNu1p375Test() {
            Complex[] expecteds = [
                "-5.27705248199702819773776225908853739551e-30+2.478879067656174354854949126828274145020e-29i",
                "-5.27705248199702819773776225908853739551e-30-2.478879067656174354854949126828274145020e-29i",
                "-9.572456252895233493029366205546736747497e26-8.19896292793433789719356411048917496803e25i",
                "-9.572456252895233493029366205546736747497e26+8.19896292793433789719356411048917496803e25i",
                "-205.2080358950215402640576240225591685438-416.7290866214760641928695100914379240115i",
                "-0.00001714958853999699091968918500505305032212-0.00004955578201524993327661980544066492052557i",
                "-205.2080358950215402640576240225591685438+416.7290866214760641928695100914379240115i",
                "-0.00001714958853999699091968918500505305032212+0.00004955578201524993327661980544066492052557i"
            ];

            foreach ((Complex z, Complex expected) in zs_mini.Zip(expecteds)) {
                Complex actual = DDoubleOptimizedBessel.ComplexBessel.BesselK(1.375, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 2e-30);
            }
        }

        [TestMethod()]
        public void HankelH1Nu1p25Test() {
            Complex[] expecteds = [
                "18.4291878726180359187566373640030564795-295.2319473632391424158850286595650387789i",
                "-2.06951583489253769330923282083229978534e-6-0.00003330852145922360428682097176358411639962i",
                "-221.7919127932079921760532896678434610892-195.7290611816043696068447244240899410570i",
                "0.00002208931271448919494119773512109508748611+0.00002501605007574010067108694777106542618294i",
                "9.11589580995419086748940818746628621207e-30+1.326373067206692081977319837757965245760e-29i",
                "-2.93296215824187966083232188070945134646e-30-1.582478564585916657645835643219594172333e-29i",
                "-5.445121799660767246044545722178884756910e26+2.819730904419201029348292537052276322562e26i",
                "5.844133392562919745278893761769145496573e26-1.856431705290732276978266742781590818347e26i"
            ];

            foreach ((Complex z, Complex expected) in zs_mini.Zip(expecteds)) {
                Complex actual = DDoubleOptimizedBessel.ComplexBessel.HankelH1(1.25, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 2e-30);
            }
        }

        [TestMethod()]
        public void HankelH2Nu1p25Test() {
            Complex[] expecteds = [
                "-2.06951583489253769330923282083229978534e-6+0.00003330852145922360428682097176358411639962i",
                "18.4291878726180359187566373640030564795+295.2319473632391424158850286595650387789i",
                "0.00002208931271448919494119773512109508748611-0.00002501605007574010067108694777106542618294i",
                "-221.7919127932079921760532896678434610892+195.7290611816043696068447244240899410570i",
                "-5.445121799660767246044545722178884756910e26-2.819730904419201029348292537052276322562e26i",
                "5.844133392562919745278893761769145496573e26+1.856431705290732276978266742781590818347e26i",
                "9.11589580995419086748940818746628621207e-30-1.326373067206692081977319837757965245760e-29i",
                "-2.93296215824187966083232188070945134646e-30+1.582478564585916657645835643219594172333e-29i"
            ];

            foreach ((Complex z, Complex expected) in zs_mini.Zip(expecteds)) {
                Complex actual = DDoubleOptimizedBessel.ComplexBessel.HankelH2(1.25, z);

                Console.WriteLine(z);
                Console.WriteLine(expected);
                Console.WriteLine(actual);

                Assert.IsTrue((actual - expected).Magnitude / expected.Magnitude < 2e-30);
            }
        }

        [TestMethod]
        public void BesselYInterpolateTest() {
            ddouble[] num16peps_expecteds = {
                "-4.971698547700163402165383494487352875503e35",
                "-7.597971175988940637927092644594564570580e30",
                "-1.158038023949481155746106008953919022306e28",
                "-1.161836491249807934322062089286378956326e26",
                "-3.273622589386874501314414461473451997375e24",
                "-1.772533393542215837165151753360079612861e23",
                "-1.506413623614907589198107199113195598227e22",
                "-1.780781482647776001596224750998972843054e21",
                "-2.708638632833416542513579804203571534922e20",
                "-5.026412180508522129740695924336459922846e19",
                "-1.095595214401546261626505441363518022984e19",
                "-2.727456011195291990849496021185648637726e18",
                "-7.591706620569295966765317045270122653990e17",
                "-2.323868238929133037407922718548632657025e17",
                "-7.721102478922891562328138202106114464966e16",
                "-2.755193376172858102691033125358059860352e16",
                "-1.046826846799931391571988645660933106627e16",
                "-4.204720609040171478858112192677769526387e15",
                "-1.774734665737716901162837586914260668779e15",
                "-7.831761067660713306196338291869143957801e14",
                "-3.597772580897664254936976501192255305850e14",
                "-1.714105604118276487267866511257309096995e14",
                "-8.442446617979542043884694856457146739572e13",
                "-4.286496132528802885036680365038854513912e13",
                "-2.238032094847674041715413070287094617767e13",
                "-1.198982166344544329466289191539002417827e13",
                "-6.578119315269153469990801523544857499605e12",
                "-3.689649635533263530733096880615949300397e12",
                "-2.112483960836015273966465194708469611309e12",
                "-1.232892720797260522787948154821513032110e12",
                "-7.325521956392611144532358047295385540820e11",
                "-4.426319698333208546124961113396423705032e11",
                "-2.717019147244114875377594393851126196177e11",
                "-1.692711802246421762220981313850733422036e11",
                "-1.069408657959333885612751359198585285770e11",
                "-6.845985173464369354500653052132980477647e10",
                "-4.437607527440752220666274429113317801045e10",
                "-2.910693520602409512746157376798747878833e10",
                "-1.930700901894948840115732947846727192325e10",
                "-1.294376941714483922273168169411675964523e10",
                "-8.766135151042003446489565820561713470574e9",
                "-5.994429825989528877732772288355373054809e9",
                "-4.136998893389975726034729705062880408514e9",
                "-2.880312709295259667367417568841143628319e9",
                "-2.022281979247201171090807870635823367174e9",
                "-1.431310590555574895078995653808138350629e9",
                "-1.020865696673530965024112057048066601375e9",
                "-7.335127524129137669607820044335568554414e8",
                "-5.307883651408199145338641320592058337885e8",
                "-3.867123051270012173520732680822916505935e8",
                "-2.835913013523384137397161071896301282574e8",
                "-2.092802588420638229963326349470571054944e8",
                "-1.553790154753685751883833221358742413631e8",
                "-1.160350569442679517562289894041092237688e8",
                "-8.714212139796302434129257352840972061483e7",
                "-6.579960741146419077245408848129493968484e7",
                "-4.994517241337513433758708140798275786788e7",
                "-3.810322987947616265566912269183009141526e7",
                "-2.921158353602887159366604699546410725065e7",
                "-2.250111866514813502639647264716904339091e7",
                "-1.741177035579478073092814644871190174632e7",
                "-1.353346158964511912269911943708236228524e7",
                "-1.056433952874654749529637075019158906239e7",
                "-8.281070323985570440275141502843528387183e6"
            };
            ddouble[] num15meps_expecteds = {
                "1.061056564180184267102667114456154744816e33",
                "3.234383574137729763154772470472785372816e28",
                "7.382954003447921522066581021384556638295e25",
                "9.865455345018435816867577371159044540268e23",
                "3.471764838724895130896834880812113385734e22",
                "2.254292159960965424820947307845894494626e21",
                "2.233939757427016605130350522247904904257e20",
                "3.016709363121320346859443112764177455978e19",
                "5.160132500736268876271068689271328090892e18",
                "1.063616117943450370585094536033813715114e18",
                "2.549472018945585012295380879269953142600e17",
                "6.922216703320809689587718153148187702239e16",
                "2.086909341352220998964574000800785187323e16",
                "6.878426444018253038032518008251373504854e15",
                "2.448283028021965543335483062372350676494e15",
                "9.317860870235719346744182526824694025071e14",
                "3.761246049050023105535904488794797050428e14",
                "1.599526338345574537625472281261864968295e14",
                "7.126099542585140437639542477966488526194e13",
                "3.310138591773604172426014111841473624521e13",
                "1.596651013789980889140091343668136253247e13",
                "7.969396930962572881377413849550583958595e12",
                "4.103714975384485376704149207584838025480e12",
                "2.174292001499835473157718412674005387814e12",
                "1.182607829179787123229084978961178674137e12",
                "6.589577969678726860268463950912201859822e11",
                "3.754748064479023692918191206586278543781e11",
                "2.184281981081631417617623328757144791409e11",
                "1.295428811708728968898203707840792294708e11",
                "7.822245989159077549960320779643855194845e10",
                "4.803458775079785176297883582901603020494e10",
                "2.996549134528829720640872879451832254341e10",
                "1.897212640572225489983723385870279708038e10",
                "1.218029552017652364965119304582098556044e10",
                "7.923186078495397377322049691976337893365e9",
                "5.218244751270392846951621863122814336665e9",
                "3.477283194264212822814400811565171624662e9",
                "2.343033565112439511946070963662190575432e9",
                "1.595484195205331416387135030572173522185e9",
                "1.097371204525293962254996352900985948415e9",
                "7.619914827644491048799547063714309657086e8",
                "5.339319733100945541378648021721210573553e8",
                "3.773797641463243007310020718826376082172e8",
                "2.689418336669507546774039048720058331716e8",
                "1.931820845223608522825620886864863024158e8",
                "1.398157764954526364281770731935786537544e8",
                "1.019267125100902070537943204418905817263e8",
                "7.482259584479950428294568675887107793150e7",
                "5.529283743894396779505546328313517921271e7",
                "4.112274972848346965487231859093489077900e7",
                "3.077267042886068694245134080025004211808e7",
                "2.316417648137810431548022998452239589594e7",
                "1.753645658428690735302625466132596173679e7",
                "1.334904321015409012273678531731551761945e7",
                "1.021543174994102469780370421914730054982e7",
                "7.857442128278242719411045519439084332074e6",
                "6.073616903818228999755573750634167298622e6",
                "4.717187985979120599414188565387707811845e6",
                "3.680612628001674632970939730926102465986e6",
                "2.884650450723245163552161403901989645146e6",
                "2.270600244264345215354160729505103836192e6",
                "1.794749487055984804123892777048962834984e6",
                "1.424382949295686953433268755250599152339e6",
                "1.134896390761379790232325322869613469374e6"
            };
            ddouble[] nup15peps_expecteds = {
                "-1.061076538646162838850178577302194003452e33",
                "-3.234444461640662447905128800445908487253e28",
                "-7.383092988086951744156978324269891700277e25",
                "-9.865641062923403912520350922156655658648e23",
                "-3.471830194947171831203444656991855553642e22",
                "-2.254334597172679667423501036255707229715e21",
                "-2.233981811503292342002924886481805094885e20",
                "-3.016766152891606969424541282666305097146e19",
                "-5.160229640600966396657181560356898738305e18",
                "-1.063636140593212249102770783056714350542e18",
                "-2.549520012939331326324078851220013098162e17",
                "-6.922347014547095362478701257283631046569e16",
                "-2.086948627570360944191526425760541415163e16",
                "-6.878555930889764271945471450413136944506e15",
                "-2.448329117125694828440526520859311049856e15",
                "-9.318036279635395783125544939146578242845e14",
                "-3.761316854776861838905083359282552399873e14",
                "-1.599556449543700717406068956907542975830e14",
                "-7.126233691920509057045061925904226876684e13",
                "-3.310200905370181207805433565267323110117e13",
                "-1.596681070859915931973711388728908042665e13",
                "-7.969546955431805809058201802706327258162e12",
                "-4.103792228113443239381394745076079041711e12",
                "-2.174332932702824632532882511243866525981e12",
                "-1.182630091857053716965961920414953298812e12",
                "-6.589702018956976405022970115450123306355e11",
                "-3.754818747880836481043932379882931703096e11",
                "-2.184323100346703240216082081046138692908e11",
                "-1.295453198249089453923328380708443328850e11",
                "-7.822393243501269192427088073969318469920e10",
                "-4.803549200535996773556600467008843083699e10",
                "-2.996605544781366078809819529938686149193e10",
                "-1.897248355736357673681740180228142625460e10",
                "-1.218052481511393793861987231110606575134e10",
                "-7.923335233042015802757487928463283996599e9",
                "-5.218342985102377708845993754949657188642e9",
                "-3.477348654369928633579275899229855333520e9",
                "-2.343077672887475190974514783495237692185e9",
                "-1.595514230309822854191403340819666604479e9",
                "-1.097391862616983532652076240087602458607e9",
                "-7.620058273088044756156281103078331829698e8",
                "-5.339420246178153582426498501934711390978e8",
                "-3.773868683474709816538438689565609841331e8",
                "-2.689468965162729732470676384411351199224e8",
                "-1.931857211889677290477305019145058100782e8",
                "-1.398184085374805124825994710914270935095e8",
                "-1.019286312877690516361073811737134808617e8",
                "-7.482400438553671287581724337553409097291e7",
                "-5.529387833057742771937992335050177416335e7",
                "-4.112352386719871483896399610466519423569e7",
                "-3.077324972659023447194196014943189116822e7",
                "-2.316461254866298692650816019964416461150e7",
                "-1.753678670934079473850827942230838063506e7",
                "-1.334929450685065571504048025941322941810e7",
                "-1.021562405617693521047363163877556709936e7",
                "-7.857590045190279518163469943231823994691e6",
                "-6.073731240092380372244103587238751071417e6",
                "-4.717276787381535274918412499485924593134e6",
                "-3.680681915798574662474255476706722445845e6",
                "-2.884704754475027342718270035711044848657e6",
                "-2.270642988476917049501345607067581851083e6",
                "-1.794783273343897774231288968455640962580e6",
                "-1.424409763406943357130261089585138738124e6",
                "-1.134917755267393763506202941843024164628e6"
            };
            ddouble[] nup16meps_expecteds = {
                "-4.971792140282170074120243448346227084166e35",
                "-7.598114208341676912679301631079516645602e30",
                "-1.158059824098400778585401467984353768625e28",
                "-1.161858362905147122754073085720252709201e26",
                "-3.273684215567085651684395611672494871126e24",
                "-1.772566761610573763949308781516937979181e23",
                "-1.506441981959495751257409745661905529250e22",
                "-1.780815005987003351381672096420896884135e21",
                "-2.708689623150093560419505442627921670972e20",
                "-5.026506803078593518140528371002504328181e19",
                "-1.095615839060094279869844011371703226460e19",
                "-2.727507355750467298055953347480760389013e18",
                "-7.591849534991326825024307657851303103674e17",
                "-2.323911985915527714563145369222705715557e17",
                "-7.721247829231243381206233056697196130156e16",
                "-2.755245242885919976194786853054134584456e16",
                "-1.046846553390458934786760164199542163173e16",
                "-4.204799763207432641051316699474829878136e15",
                "-1.774768075245181905696776337559150186400e15",
                "-7.831908501124940691601952180085785547063e14",
                "-3.597840309224731628025585110473365143963e14",
                "-1.714137872279304823803420733443757367123e14",
                "-8.442605547643086450220937512052837316793e13",
                "-4.286576826126163200908708329463931722178e13",
                "-2.238074225962495452084137924870073259224e13",
                "-1.199004737269883326216291575853844821540e13",
                "-6.578243148837411303721622830480895997452e12",
                "-3.689719093452478941756190006885491642743e12",
                "-2.112523728498202633067434790959958469470e12",
                "-1.232915930091221039380229281041582477627e12",
                "-7.325659859869233335734409456545365506098e11",
                "-4.426403024119231038172802017685782998061e11",
                "-2.717070295324577190405111788209113041317e11",
                "-1.692743667667594309975619566037089148428e11",
                "-1.069428789654077325671882595038347821885e11",
                "-6.846114049626633922699684807200332399308e10",
                "-4.437691065720967851954624032790698224987e10",
                "-2.910748314616851375275665216615686249506e10",
                "-1.930737247478004061234488570327547913576e10",
                "-1.294401308453316267126285153761277609776e10",
                "-8.766300174165415664346112274288111106354e9",
                "-5.994542671561256409695364958549566999664e9",
                "-4.137076772691079012374914503842004992815e9",
                "-2.880366931388748074468007017990546331754e9",
                "-2.022320048850610827216245547445186275908e9",
                "-1.431337535080175156434731191473569417509e9",
                "-1.020884914543544003868573647574504538920e9",
                "-7.335265608431261945477050293273179180002e8",
                "-5.307983572698912979523890033538555380709e8",
                "-3.867195850139017074386281017441203812116e8",
                "-2.835966399789412827971934881246815396520e8",
                "-2.092841985579577228973390330444183737185e8",
                "-1.553819404964873076285784331404021074526e8",
                "-1.160372413125434033990318258695180597980e8",
                "-8.714376185465301150724598862010245320662e7",
                "-6.580084609379651219599088301911375943888e7",
                "-4.994611263483065652031061787511000798653e7",
                "-3.810394717551352268290196457673628231443e7",
                "-2.921213344618561374220878013775468247924e7",
                "-2.250154225032203578165268584279850269660e7",
                "-1.741209813362146731921290584001244299192e7",
                "-1.353371635803091533820628844516595655171e7",
                "-1.056453840319643703191195946542991494290e7",
                "-8.281226215728797528494400865472624544762e6"
            };
            ddouble[] nu0peps_expecteds = {
                "-1.839534518413228609419738569348901487959",
                "-1.392731165574072441778454060319215243090",
                "-1.126956113001005191407251658399948150301",
                "-0.9345876595295288738098174581162201008103",
                "-0.7821176846434801920018733980334970128311",
                "-0.6547565189231948188720339330635825515956",
                "-0.5447094276782058614268776417118601970051",
                "-0.4473941084068852779530530736694012215845",
                "-0.3599097366178956527248614188276403213752",
                "-0.2803196114230786865169882096154731035453",
                "-0.2072782883098642905205034539335428514162",
                "-0.1398220760108598952820523544371704245376",
                "-0.07724393679303663902894032009964915744269",
                "-0.01901503338866927823176923986340620094711",
                "0.03526650009518621806397063557582652981537",
                "0.08591016568988390438550310745832944093803",
                "0.1331571674694380539392455118969573755934",
                "0.1771975705633272849041089061229348870142",
                "0.2181828050085699534561338171116989011792",
                "0.2562349531486516884494564506498099210159",
                "0.2914537650010177453174249016012072716985",
                "0.3239220365188048354762975421382796856698",
                "0.3537097867623044195842434923389135950403",
                "0.3808775391531367015570853607948217292446",
                "0.4054789241050808145859685362072249533475",
                "0.4275627601843576408165952886742166320520",
                "0.4471747290887033463302447306785731759901",
                "0.4643587301407676812452770756001191476226",
                "0.4791579787748234939504962643583066214478",
                "0.4916158980883831342918976885481013962602",
                "0.5017768412091533620762293185514638528011",
                "0.5096866738189686641993942830599912624039",
                "0.5153932398686869909804520056492576609858",
                "0.5189467287434180762696641192244871048741",
                "0.5203999584931350960472050236885779248951",
                "0.5198085869408739676558457280038348364471",
                "0.5172312603099932522671371328045219918196",
                "0.5127297073200578993251337267535567639276",
                "0.5063687853745108950120789251930022470070",
                "0.4982164844179000613687446674567426361782",
                "0.4883438932124721873359743999206166388586",
                "0.4768251321252448337095292056966432291379",
                "0.4637372559904059450422602124996322566537",
                "0.4491601301897495469127146040386820382892",
                "0.4331762827538909704542634051395246820990",
                "0.4158707350121094244799217349288904489292",
                "0.3973308130954180807336541295960933287444",
                "0.3776459424152716362104535265879356874678",
                "0.3569074270907857832093428442512106220359",
                "0.3352082161737607876079349184133470343273",
                "0.3126426584178023813811258402865004077071",
                "0.2893062472510952509862698083270807474996",
                "0.2652953575384150366362125662788618225152",
                "0.2407069756539320774726755260899681759957",
                "0.2156384243299613284534814937941043298070",
                "0.1901870836961711978602646710543634703032",
                "0.1644501098773478524008550937004863165511",
                "0.1385241524743702141535441417697526120543",
                "0.1125050722115640499684393156327690238494",
                "0.08648765999324048890303079747782805445322",
                "0.06056535857231136095663797729714438638221",
                "0.03482998799386818678569395989816495393367",
                "0.009371475936080388813019609031574708487798",
                "-0.01572240597062933387379850363029957758946"
            };
            ddouble[] nu0meps_expecteds = {
                "-1.833404488447961204612316434334274283917",
                "-1.386619131581441106674082534244301310877",
                "-1.120873963886776151805272279406182005928",
                "-0.9285472112295295655274129803004905069106",
                "-0.7761306350629374805903211976107499646490",
                "-0.6488344113817084460402185846091374488793",
                "-0.5388636164168292104032686860151133140324",
                "-0.4416357251793123411123250623208094299342",
                "-0.3542496581843939983846383118382691769977",
                "-0.2747684279531336907153525786253948905487",
                "-0.2018462727555452834638298081523687104160",
                "-0.1345191545315041308244780840401972089227",
                "-0.07207966031745942839720365741483553539642",
                "-0.01399855040148357209873447727546213847871",
                "0.04012646945860777700195813499272863023864",
                "0.09060535415860412456621995088758709124719",
                "0.1376797836852116419950393187064682597184",
                "0.1815403205894494411234195394706999253577",
                "0.2223389122318633936945306502882946486380",
                "0.2601981765124633386689306469993007747655",
                "0.2952184155126784813878651747724510633162",
                "0.3274829919808062143013209861474255286807",
                "0.3570625046824338119639592494053518888666",
                "0.3840180677949588435761033640795204332099",
                "0.4084039116464338313585163844271036049219",
                "0.4302694619579204397330909937921423985951",
                "0.4496610128851648783198382684291648973837",
                "0.4666230795616631887706768207996765372041",
                "0.4811994946286749447129141052327305399254",
                "0.4934342978303671540013644712065974911104",
                "0.5033724564310968355904466848386738179544",
                "0.5110604458030166424537671563608906121157",
                "0.5165467132244767027466605606543186997578",
                "0.5198820431550040493868658239574758952470",
                "0.5211198386082213269841880843796331994238",
                "0.5203163304410978284181304162426175334272",
                "0.5175307242070509501627150883751684806737",
                "0.5128252925283419121960667439487215820254",
                "0.5062654196166256105769569838033894753216",
                "0.4979196035249132557662287928411841452074",
                "0.4878594208860281366006700819290024967869",
                "0.4761594582336917298221559526732669395558",
                "0.4628972134758587300669462435245723117353",
                "0.4481529706675030662504545618757413172048",
                "0.4320096508898022542262517865105081230466",
                "0.4145526417674677551270174655507364140750",
                "0.3958696079324042104452493114176708588220",
                "0.3760502845593611349357543147524374634165",
                "0.3551862559493652626666396912980858893638",
                "0.3333707210127901122731402601647271424981",
                "0.3106982474005626470939004627341734211492",
                "0.2872645159449042663936928815535067496930",
                "0.2631660569966635639132603620570711196558",
                "0.2384999801818907777576234169622572746519",
                "0.2133636990435271485830274312702313677513",
                "0.1878546519830611048774294640023826610445",
                "0.1620700208702069939287927115481495071567",
                "0.1361064486448417247988445446758903789565",
                "0.1100597571935711453843125757395944884079",
                "0.08402466674256182887022016066109590478929",
                "0.05809451796799287453529522329846446243278",
                "0.03236099798511488502295945175261277070959",
                "0.006913871336017663075667988201013925585490",
                "-0.01815928294553505601749542527148216115359"
            };
            ddouble[] nu4peps_expecteds = {
                "-2.021845854676099261585501971434324049030e6",
                "-126317.6746765963772357557412470793021700",
                "-24972.51639920326990273597351833480401929",
                "-7915.062074558463435379279198330235704604",
                "-3250.130589415402226674267760036198183407",
                "-1572.472730255684445689618868987682765715",
                "-852.1477771017682640076062967987844768611",
                "-501.8492849647438365082342920620133716164",
                "-314.9886191503584936449729112036802558722",
                "-207.9227243106474907036606101918246533769",
                "-142.9798698937967810235888729037753240663",
                "-101.7122315545516126505439825783403550989",
                "-74.45420590832517322900708656646270716613",
                "-55.85098525139777156724401352304439463347",
                "-42.79359905813413276845002189611557824421",
                "-33.40340437870890083073192075470675858751",
                "-26.50536925688959512576821858505550257512",
                "-21.34210680467178337667257758645404061626",
                "-17.41233479194516341564499401761583071315",
                "-14.37639230577409506027501690583682764333",
                "-11.99925105998912826527698767215629663758",
                "-10.11518697814326179963934461776035915068",
                "-8.605334485353697077454754824775467756081",
                "-7.383106070724945982330774966715419237396",
                "-6.384528699234892078843304783714736120660",
                "-5.561720615046652732749340324142225800858",
                "-4.878413585463394389564151848198784906298",
                "-4.306831529839975377940071932980268896367",
                "-3.825483568755183027713826649836682412473",
                "-3.417582998696550371057008689722732731766",
                "-3.069900799776075532014828083500048044280",
                "-2.771924782491753003684033978069461645601",
                "-2.515236349608786332843010098258942945746",
                "-2.293043973307885396803783207212106790085",
                "-2.099830738731040762269368276877227096312",
                "-1.931085744751218673398347266358978829404",
                "-1.783097734488453585866862997543566648867",
                "-1.652795315747727145606333961245092954803",
                "-1.537622354049421225366644260962345636847",
                "-1.435440128696376597771097255537487115444",
                "-1.344450005244519639870973107314687203652",
                "-1.263131947144305458930247881871592689398",
                "-1.190195337784954557232565701424710159073",
                "-1.124539431407426087281226221742843326581",
                "-1.065221381170990410689823896002031897744",
                "-1.011430264278625132576550900031537196205",
                "-0.9624658796842490865748804238372348122170",
                "-0.9177213638333101502178368290958880341324",
                "-0.8766688760738250169526488911424762956029",
                "-0.8388477638245007149074826068302865812909",
                "-0.8038547400576970850282794460399015797908",
                "-0.7713357008427422247044349780654351350715",
                "-0.7409788850713940110690174677663407071756",
                "-0.7125091368992031012501380853889905477143",
                "-0.6856830775361220003152646806076605452702",
                "-0.6602850295741253990134990764047571466481",
                "-0.6361235661571956526112634977603697842349",
                "-0.6130285805946200241229924165021875739070",
                "-0.5908487907352825600503013693894161300982",
                "-0.5694496075201964743736657353676455484002",
                "-0.5487113093603716347894138748539292060226",
                "-0.5285274739308416536650647461368610584684",
                "-0.5088036270869385175394779416222995572074",
                "-0.4894560752561055353643654747398985825527"
            };
            ddouble[] nu4meps_expecteds = {
                "-1.984896094179859796758692895228798169621e6",
                "-124345.5675676393049860264519250935695910",
                "-24621.65683057359871890462636522350386104",
                "-7812.654735908085986386339260674892490680",
                "-3210.889421591162779255159335814283286604",
                "-1554.601291888909576651868905180099813292",
                "-842.9751701967680107141026704212995681906",
                "-496.7096137378086264637830003675831384940",
                "-311.9084971786970250634213534716205817812",
                "-205.9760585513782701214752938627440156459",
                "-141.6953257488269590773096191021405182867",
                "-100.8337774708640933323037851517546239444",
                "-73.83512010882928589829815070536482836626",
                "-55.40333787905857805783366747163512223994",
                "-42.46265176042458599226665543851652672725",
                "-33.15394807949779846404231540307947893709",
                "-26.31410292178795098137731412694909470549",
                "-21.19321990793064029192863077186663346715",
                "-17.29485943138058300791960805255163268311",
                "-14.28256837864515334396352895970001708188",
                "-11.92348999539812950110677584548070021023",
                "-10.05339921753872390299006078709767412075",
                "-8.554483668575686571771795744729286973512",
                "-7.340907771637680242926484441343630130548",
                "-6.349243251522629771110713616180110981749",
                "-5.532008347096776224782774185074051212100",
                "-4.853232169976874819392126537514851961914",
                "-4.285362206627699091183693513865407852668",
                "-3.807077453560446770267597358134693319716",
                "-3.401721555076528764547685522554816082397",
                "-3.056166423541574008381810178241403278971",
                "-2.759978691275948968392870509655771080360",
                "-2.504801807236456597040003646678992294351",
                "-2.283893449901105523348104990608898621084",
                "-2.091776002212387230072181764809376925350",
                "-1.923970155340596332595882341504861436604",
                "-1.776790208802334739100489792487854265684",
                "-1.647185566354484635606320863827726636913",
                "-1.532617110726210224420701447393799322103",
                "-1.430960120676772124101426378271371175991",
                "-1.340427537352615061907922980924519848189",
                "-1.259508942375473077456326034289399180027",
                "-1.186921748470387309994125694404213185825",
                "-1.121571943317363848384132208191374726561",
                "-1.062522351718467345054539270107778743909",
                "-1.008966848784413435766293994593721588279",
                "-0.9602093094723861996848178097512249662807",
                "-0.9156463474894746051285665235143995174867",
                "-0.8747531010639758061793926640530567170137",
                "-0.8370714802450164056925791934595961533350",
                "-0.8022004118740507691839397307477510942452",
                "-0.7697877127993135102542839102836587513667",
                "-0.7395232956924481840805168229688679595221",
                "-0.7111334697811344370345492978959107050727",
                "-0.6843761445536850078430450898765266257149",
                "-0.6590367807653184674115367281453359167505",
                "-0.6349249619719746007077275520693713923586",
                "-0.6118714829376798911357708406218055339720",
                "-0.5897258698385929230441939721382118829401",
                "-0.5683542621748719275246455669880213037430",
                "-0.5476375984418760797303213563531868831367",
                "-0.5274700574838946094774805813756979902334",
                "-0.5077577155106374466260023114960020346631",
                "-0.4884173853567055126315035181794709977541"
            };

            for (int n = -16; n <= 16; n++) {
                foreach (ddouble u in new ddouble[] {
                    Math.ScaleB(-1, -8), Math.ScaleB(-1, -16), Math.ScaleB(-1, -20), Math.ScaleB(-1, -25), Math.ScaleB(-1, 96),
                    Math.ScaleB(1, -8), Math.ScaleB(1, -16), Math.ScaleB(1, -20), Math.ScaleB(1, -25), Math.ScaleB(1, -96) }) {
                    ddouble nu = n + u;

                    if (ddouble.Abs(nu) > 16) {
                        continue;
                    }

                    for (ddouble x = 1d / 32; x <= 8; x += 1d / 32) {
                        Complex y = DDoubleOptimizedBessel.ComplexBessel.BesselY(nu, x);
                        Complex y_dec = DDoubleOptimizedBessel.ComplexBessel.BesselY(ddouble.BitDecrement(nu), x);
                        Complex y_inc = DDoubleOptimizedBessel.ComplexBessel.BesselY(ddouble.BitIncrement(nu), x);

                        Console.WriteLine($"{nu}, {x}");
                        Console.WriteLine(y);
                        Console.WriteLine(y_dec);
                        Console.WriteLine(y_inc);

                        Assert.IsTrue(ddouble.Abs((y_dec.R - y.R) / y.R) < 8e-28);
                        Assert.IsTrue(ddouble.Abs((y_inc.R - y.R) / y.R) < 8e-28);
                    }

                    for (ddouble x = 1d / 32; x <= 8; x += 1d / 32) {
                        Complex y = DDoubleOptimizedBessel.ComplexBessel.BesselY(nu, (x, x));
                        Complex y_dec = DDoubleOptimizedBessel.ComplexBessel.BesselY(ddouble.BitDecrement(nu), (x, x));
                        Complex y_inc = DDoubleOptimizedBessel.ComplexBessel.BesselY(ddouble.BitIncrement(nu), (x, x));

                        Console.WriteLine($"{nu}, {(x, x)}");
                        Console.WriteLine(y);
                        Console.WriteLine(y_dec);
                        Console.WriteLine(y_inc);

                        Assert.IsTrue((y_dec - y).Magnitude / y.Magnitude < 8e-28);
                        Assert.IsTrue((y_inc - y).Magnitude / y.Magnitude < 8e-28);
                    }

                    for (ddouble x = 1d / 32; x <= 8; x += 1d / 32) {
                        Complex y = DDoubleOptimizedBessel.ComplexBessel.BesselY(nu, (0, x));
                        Complex y_dec = DDoubleOptimizedBessel.ComplexBessel.BesselY(ddouble.BitDecrement(nu), (0, x));
                        Complex y_inc = DDoubleOptimizedBessel.ComplexBessel.BesselY(ddouble.BitIncrement(nu), (0, x));

                        Console.WriteLine($"{nu}, {(0, x)}");
                        Console.WriteLine(y);
                        Console.WriteLine(y_dec);
                        Console.WriteLine(y_inc);

                        Assert.IsTrue((y_dec - y).Magnitude / y.Magnitude < 8e-28);
                        Assert.IsTrue((y_inc - y).Magnitude / y.Magnitude < 8e-28);
                    }

                    for (ddouble x = 1d / 32; x <= 8; x += 1d / 32) {
                        Complex y = DDoubleOptimizedBessel.ComplexBessel.BesselY(nu, (-x, x));
                        Complex y_dec = DDoubleOptimizedBessel.ComplexBessel.BesselY(ddouble.BitDecrement(nu), (-x, x));
                        Complex y_inc = DDoubleOptimizedBessel.ComplexBessel.BesselY(ddouble.BitIncrement(nu), (-x, x));

                        Console.WriteLine($"{nu}, {(-x, x)}");
                        Console.WriteLine(y);
                        Console.WriteLine(y_dec);
                        Console.WriteLine(y_inc);

                        Assert.IsTrue((y_dec - y).Magnitude / y.Magnitude < 8e-28);
                        Assert.IsTrue((y_inc - y).Magnitude / y.Magnitude < 8e-28);
                    }

                    for (ddouble x = 1d / 32; x <= 8; x += 1d / 32) {
                        Complex y = DDoubleOptimizedBessel.ComplexBessel.BesselY(nu, -x);
                        Complex y_dec = DDoubleOptimizedBessel.ComplexBessel.BesselY(ddouble.BitDecrement(nu), -x);
                        Complex y_inc = DDoubleOptimizedBessel.ComplexBessel.BesselY(ddouble.BitIncrement(nu), -x);

                        Console.WriteLine($"{nu}, {-x}");
                        Console.WriteLine(y);
                        Console.WriteLine(y_dec);
                        Console.WriteLine(y_inc);

                        Assert.IsTrue((y_dec - y).Magnitude / y.Magnitude < 8e-28);
                        Assert.IsTrue((y_inc - y).Magnitude / y.Magnitude < 8e-28);
                    }
                }
            }

            ddouble eps = Math.ScaleB(1, -9);

            foreach ((ddouble nu, ddouble[] expecteds) in new (ddouble, ddouble[])[] {
                (-16 + eps, num16peps_expecteds), (-15 - eps, num15meps_expecteds), (15 + eps, nup15peps_expecteds), (16 - eps, nup16meps_expecteds),
                (+eps, nu0peps_expecteds), (-eps, nu0meps_expecteds), (4 + eps, nu4peps_expecteds), (4 - eps, nu4meps_expecteds),
            }) {
                for ((ddouble x, int i) = (1d / 16, 0); i < expecteds.Length; x += 1d / 16, i++) {
                    Complex expected = expecteds[i];

                    Complex actual = DDoubleOptimizedBessel.ComplexBessel.BesselY(nu, x);

                    Console.WriteLine($"{nu}, {x}");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    Assert.IsTrue(ddouble.Abs((actual.R - expected.R) / expected.R) < 8e-28);
                }
            }
        }

        [TestMethod]
        public void BesselKInterpolateTest() {
            ddouble[] nu15peps_expecteds = {
                "1.666502653570570555725780557421577778568e33",
                "5.077819479651207054123048417111991584555e28",
                "1.158278512587413236267055560315869340522e26",
                "1.546236481189214991378981117639017751544e24",
                "5.434553460125745649824705504637834044988e22",
                "3.523363012692633368790800559894090575298e21",
                "3.485227353894452177986526187256486535788e20",
                "4.696609314419667625598225302206374624888e19",
                "8.014602091023147975163582890881239434048e18",
                "1.647612243794438080889201313930376919441e18",
                "3.937750360278051856110696957412876848416e17",
                "1.065736379359135204865214139143282666418e17",
                "3.201796259661618613367410778058047597299e16",
                "1.051340893403686855740538856411928364713e16",
                "3.726998957313661758067204914902866571808e15",
                "1.412329073212969512516651458776528190203e15",
                "5.674822596388879852508331448579103702163e14",
                "2.401550068293613791573005096616869930561e14",
                "1.064414240235128885340372325695523199087e14",
                "4.917477109571985725876916189320454650040e13",
                "2.358426702318160176078033802680251412966e13",
                "1.170126693858276226515595828307009509835e13",
                "5.987678816119117723273902972071503524060e12",
                "3.151750892498314659375435285415493041987e12",
                "1.702575201497018103435207713016797736274e12",
                "9.419622932087969512641087981483752655187e11",
                "5.327775380781879313557583814016673044066e11",
                "3.075684725945137024534584332325698209338e11",
                "1.809645431240555674335468440867731990801e11",
                "1.083769807987248964431547109059214435412e11",
                "6.598786802925902799387129393765445898156e10",
                "4.080512280121877444786519483605945533198e10",
                "2.560185027513340809325589329866130507422e10",
                "1.628373489151774437154400221105351748979e10",
                "1.049097127505303033611993503279184471977e10",
                "6.841301271545929878053624149442710722554e9",
                "4.512648832805159682663988550144395512941e9",
                "3.009026160048775522614046278347870759321e9",
                "2.027097780266034824696441296888837021385e9",
                "1.378952299295672977415332250998654469443e9",
                "9.467559348855769164499254280070198313132e8",
                "6.557598421418347623973201779450126026745e8",
                "4.580228987324514585986295050718165510198e8",
                "3.224744780915794540952051215043575163536e8",
                "2.287762509131234038629352848070220237846e8",
                "1.634880670554845478879563832695042783062e8",
                "1.176474065579543173949756148459003242638e8",
                "8.522563176060218352240268743932809783604e7",
                "6.213382528221253966910204294110400319089e7",
                "4.557658079906002828094387969374135512485e7",
                "3.362821096315187990517583038131057386726e7",
                "2.495243945939370087518871885734642761743e7",
                "1.861547496216217143965028839310555821581e7",
                "1.396037515331405427205550220453277286865e7",
                "1.052196537551257136772609351955881930392e7",
                "7.968801665454103049424845234086542381881e6",
                "6.063304436015830017695987479003897776520e6",
                "4.634189146349801971779821740268172105187e6",
                "3.557277023084080624472575866099512601979e6",
                "2.742054951050554767752565740558826030918e6",
                "2.122205183665222855160357571937829413016e6",
                "1.648893380047760400952309905839121627372e6",
                "1.285984775052639223656146799267833895645e6",
                "1.006616026443457157523056861510975536509e6"
            };

            ddouble[] nu16meps_expecteds = {
                "7.808655880853390300248736248115437219773e35",
                "1.192887450581434876092553807227903072372e31",
                "1.816945359366853822418413468145286467048e28",
                "1.821244139829754231091283800959385477391e26",
                "5.125576954277155531814799204421131149318e24",
                "2.771318608006693788118694040734928704554e23",
                "2.351262064901317608272614984990368121262e22",
                "2.774080703404504867095225666545095979729e21",
                "4.210154897626975230034414136381387729966e20",
                "7.793464084901443480217591068269507257542e19",
                "1.694083788483651278224506657165284563422e19",
                "4.204764945829420589376389167397062997436e18",
                "1.166566443755447377767140534219562263243e18",
                "3.558398199789522780696245221412723903532e17",
                "1.177828967659151478457266276713540681343e17",
                "4.186024400992544679005959585332325274155e16",
                "1.583646201020631467517052520142619226843e16",
                "6.332001102641638601483453954441629011546e15",
                "2.659773670969749605820843884789523386858e15",
                "1.167790602189393782752832901086240598684e15",
                "5.336056677834474368502036430282417215897e14",
                "2.528089025688243517224574784917230000305e14",
                "1.237878601993842281621107235501139338989e14",
                "6.246748095152682240289697867166016546781e13",
                "3.240757765903198975195114363905932894438e13",
                "1.724680147426785759889644222445993416020e13",
                "9.397234053023221892902084818609046260482e12",
                "5.233265343519821337385759609027079797679e12",
                "2.974111520440662436563777065929046689012e12",
                "1.722472068767849260535705213426428236083e12",
                "1.015349153912857179825705552671510003660e12",
                "6.084941617800832439714813830823607733932e11",
                "3.703650497733049138591096744212685687187e11",
                "2.287340592115151266706111076204635691552e11",
                "1.432151229793406458897929170046909268363e11",
                "9.083753829622528041426176274739237099000e10",
                "5.832427549216462499744145465327649530816e10",
                "3.788391395702616038089397245718075000610e10",
                "2.487815068283298284766636570151745038396e10",
                "1.650803651855706635945380163356923523890e10",
                "1.106270149978566457254397900631449233444e10",
                "7.483523118098374819319695429784774571035e9",
                "5.107820496425685890756824494778111739595e9",
                "3.516160231014865912856055744543360684253e9",
                "2.440261187990150456319889454630910875613e9",
                "1.706792983300142149616806177795151708281e9",
                "1.202692790667182386706068884491329291617e9",
                "8.535317034621719226594531463851226981429e8",
                "6.098822898175068016214620069024275790672e8",
                "4.386437983968570309588446488486251040776e8",
                "3.174702286324909561865339291047856169996e8",
                "2.311591926721321342652651312418639243912e8",
                "1.692913510535518270585057052614453151484e8",
                "1.246745487079970903132040679486293984333e8",
                "9.231017104437618213288230743155689020118e7",
                "6.870118337423672210662732666718742077593e7",
                "5.138549535631185779162070190002387350538e7",
                "3.861904636876048913256361727115436212332e7",
                "2.915910345234612777922727020846302472768e7",
                "2.211510750925317657828802231849772743392e7",
                "1.684536007108334136928664605270630588632e7",
                "1.288502324247313232195966491999466314668e7",
                "9.895644886672748678585166180261674222757e6",
                "7.629571258744425015051472670610558095755e6"
            };

            ddouble[] nu0peps_expecteds = {
                "2.892341579688986530034369150974401303018",
                "2.207881316562347546988298747898686650313",
                "1.814500204058855423369075478013720263195",
                "1.541512442193464917570011970496527355084",
                "1.335143637065933360606013948149450199041",
                "1.171284261942179593344991435001473456261",
                "1.036976402556370496162654828400469906118",
                "0.9244212277142615914713743668934173913319",
                "0.8285442323513232073812295242815106780101",
                "0.7458493739582825342125154581570010919569",
                "0.6738205472638626317742370731530001916741",
                "0.6105834833324062851823229880037508132935",
                "0.5547026969084724933720212895833624608968",
                "0.5050533725834092650557288102567815064796",
                "0.4607371785261711492712485617922004384097",
                "0.4210250253440369410373118243574384254159",
                "0.3853169975339611243470470197159026835306",
                "0.3531135932288979927018077681610173974153",
                "0.3239946244230232416406414748069846977799",
                "0.2976034367067040456757557485602397276788",
                "0.2736349045599352929423571410539393908744",
                "0.2518261589896272533047562547483822019529",
                "0.2319493272778887775059712457202403018414",
                "0.2138057778865869690501531715523482755571",
                "0.1972215074031273160095066250803500894382",
                "0.1820434052887746989079901891141125755827",
                "0.1681362013468228996864023866036700874976",
                "0.1553799499677775961820879629363393745386",
                "0.1436679406355642347847713173644336780103",
                "0.1329049500622978488046813678329711753239",
                "0.1230057704645147204703788365787997301602",
                "0.1138939628176214443280433664242325189538",
                "0.1055007947561547891821939913495256496423",
                "0.09776433105841095403789298570997219871559",
                "0.09062865102870930420300158083272515157768",
                "0.08404317204641867486669085632403103049033",
                "0.07796206243508342897032657575043750790615",
                "0.07234372987259628074111059134459385578330",
                "0.06715037400353554404295894498871125606715",
                "0.06234759386906708984869791552787530509909",
                "0.05790404234502866217035997494558713369422",
                "0.05379112105622660837582365625464770961944",
                "0.04998271027686584293883166980639976920776",
                "0.04645492918148392873136215949801983204662",
                "0.04318592251513968419690117170802296527454",
                "0.04015567033514929592926504904786994959324",
                "0.03734581796237818178731946701482432318955",
                "0.03473952368620209071542045483178767564021",
                "0.03232132210825010788379252579883946358466",
                "0.03007700129754480208539940824487221593176",
                "0.02799349217300091590243574842393451960469",
                "0.02605876873598516568289615625109078062083",
                "0.02426175795190979185110201827803052704384",
                "0.02259225823063644348120448177103627247623",
                "0.02104086558490628465499973739204737358106",
                "0.01959890665746196962081489438779340237850",
                "0.01825837790377581266610078773927602627080",
                "0.01701189030065874349197701454380077465137",
                "0.01585261902342165738315323966943791902048",
                "0.01477425759730286398411102664133988112764",
                "0.01377097608391076154717011027280164749955",
                "0.01283738291159360795843791199890159230260",
                "0.01196849000089516772520901541439651995558",
                "0.01115968087339595180095175559206121583083"
            };

            ddouble[] nu4meps_expecteds = {
                "3.115836971180018751854559639041288969013e6",
                "194813.2418265801268750907894781549533872",
                "38449.50841713358174784369314237549165243",
                "12144.83502974170503620903754015084781510",
                "4962.172638335667508405772440510865365786",
                "2385.355021301430641684442262979623500716",
                "1282.537317646750795815928282624692641425",
                "748.3632157310226198150828360495397843325",
                "464.7555876401261168898926963981183787059",
                "303.1334675188553680917345948611821070067",
                "205.6949289583500637020459640156682818597",
                "144.1965540560733944406543582854577503029",
                "103.8773313907497296787562050230923043412",
                "76.58235747091924224380271526323537791420",
                "57.59149810859731392594966834320279763376",
                "44.06211609802979064350253824305315971025",
                "34.22269907394529594357207795845015704692",
                "26.93585981490422934528152098419787968708",
                "21.45196618663154433542615419389866941376",
                "17.26534923633274976219740370903001410759",
                "14.02777179704102732904086330153613420567",
                "11.49491674523927027023461728794884665625",
                "9.492476615116394741577072896074401242560",
                "7.894186150842356760663050253068900525023",
                "6.607304779305543888360916461381008647738",
                "5.562846580028823930742901108549565866706",
                "4.708895099783916168477677188051052929633",
                "4.005958703832216742672090848883951850303",
                "3.423697983583768289247379956109004294191",
                "2.938589783997985889300394506806968228538",
                "2.532239608407360288120455756276026325618",
                "2.190148726861679020889171536345408592510",
                "1.900804039418824964536381024426327607757",
                "1.654999630589657260653121226039908520467",
                "1.445326405346598016052099648835989606654",
                "1.265784870163888432190127897347715794233",
                "1.111488976293682323108025679836002162491",
                "0.9784378905859906874921777978438529151198",
                "0.8633388546486494053946796431322787220601",
                "0.7634687667618681571071238563313741119670",
                "0.6765653300111856100496870142570438803339",
                "0.6007409326176287130809072405602399853736",
                "0.5344141215618571947159068379105139901486",
                "0.4762547777274399204489230958513841418313",
                "0.4251400253100647141570669801563818179019",
                "0.3801185985675529856247943861372575719915",
                "0.3403819079958476673894421632925954094720",
                "0.3052404407862880435302763429784276908336",
                "0.2741044295169101383595871737130647694612",
                "0.2464679521578886101938062671797813771938",
                "0.2218958030009274400998272581222611448076",
                "0.2000126108728224900335653087738190704252",
                "0.1804937874837024989454380123294957564528",
                "0.1630579721032738439381562662905024793259",
                "0.1474607042997359898545170593536789920410",
                "0.1334891082568229238984380591917291631820",
                "0.1209574132757567687992118370298130647240",
                "0.1097031678170232916765781237612532899128",
                "0.09958403064330680526231750328832381478206",
                "0.09047504367960675431072894879515996373732",
                "0.08226630818692270819453383425085648328772",
                "0.07486099959086180771068693619956139813378",
                "0.06817366747214664047730756200419207079435",
                "0.06212877632747750097547154578943554977045"
            };

            ddouble[] nu4peps_expecteds = {
                "3.173842410624172802001445648783485859193e6",
                "197903.6296081332218424693477645357838771",
                "38997.71294362382092995699303766696468110",
                "12304.19496312844035967100346555604471858",
                "5022.923188693009198501793216938460334097",
                "2412.850324605925241236099047622766470471",
                "1296.546789310211542150817787671620714047",
                "756.1479218969387138130187044272427907651",
                "469.3773687164490676586547134758005552481",
                "306.0243664496637366323955797742419229594",
                "207.5810308113845723482739380855829424913",
                "145.4706249440939143613203754364092754932",
                "104.7634202133608318853991901573008202017",
                "77.21406762887029320530260244811461649185",
                "58.05155345711010585200956050435736680847",
                "44.40341915274588772966806347279800009874",
                "34.48004380652781324866775026112878785544",
                "27.13269965155442485638470275809845574163",
                "21.60445729198786807739269394000480154762",
                "17.38483805240055687554033047043699864124",
                "14.12236592396990436944776049157912086149",
                "11.57050038769012695227896118712464327784",
                "9.553380921392587044293266596759830928853",
                "7.943640221838221439388838080742113728317",
                "6.647744397162759062272862231485237209088",
                "5.596128900325740967089601251015218596574",
                "4.736450123077022891183552060778934298368",
                "4.028897573010232466027313668089599802174",
                "3.442891368528898933559380716238358856603",
                "2.954725263433533548679242224366481245689",
                "2.545864121511015724590665180926318788196",
                "2.201700288838021874439399789841627991574",
                "1.910635653924124159788829693887244174414",
                "1.663397473710077367196503199105131038714",
                "1.452523764674253802739738841653880262916",
                "1.271972910579296579805558789069625840658",
                "1.116825113891955667396403209902908797892",
                "0.9830523496967588908914951977366353664805",
                "0.8673398345492611002494981427082047142260",
                "0.7669465371350258826574933366075768050940",
                "0.6795954955620283451998492730878482326529",
                "0.6033870493031864682383509276927334123600",
                "0.5367298038345470227816320474040811111340",
                "0.4782854039852359197316434564617481793139",
                "0.4269241244819609812589967970595927097480",
                "0.3816889823639160781815962012435542101578",
                "0.3417665993141774315581693775127898819860",
                "0.3064634379866473337812785726879491418748",
                "0.2751863379623944175578509360651914801028",
                "0.2474265079586552216129487372583661425197",
                "0.2227463088626143559541859320533550159991",
                "0.2007683000001508430023752663437541433574",
                "0.1811661283780609969777408257425707932139",
                "0.1636569246301043578605204094493248621375",
                "0.1479949354436602436437311579011686748266",
                "0.1339661744195202211527222385246716290026",
                "0.1213839147188929384132903040475935633721",
                "0.1100848798445516161453723507357003001283",
                "0.09992601530335794486346449372783287865537",
                "0.09078174510628983543775506876858358273199",
                "0.08254163416560020495378316933231387170848",
                "0.07510839149239341533007358618202519692347",
                "0.06839616034291734164660478246901239528474",
                "0.06232905062742706209315361647640592093629"
            };

            for (int n = -16; n <= 16; n++) {
                foreach (ddouble u in new ddouble[] {
                    Math.ScaleB(-1, -8), Math.ScaleB(-1, -20), Math.ScaleB(-1, -25), Math.ScaleB(-1, 96),
                    Math.ScaleB(1, -8), Math.ScaleB(1, -20), Math.ScaleB(1, -25), Math.ScaleB(1, -96) }) {
                    ddouble nu = n + u;

                    if (ddouble.Abs(nu) > 16) {
                        continue;
                    }

                    for (ddouble x = 1d / 4; x <= 8; x += 1d / 4) {
                        Complex y = DDoubleOptimizedBessel.ComplexBessel.BesselK(nu, x);
                        Complex y_dec = DDoubleOptimizedBessel.ComplexBessel.BesselK(ddouble.BitDecrement(nu), x);
                        Complex y_inc = DDoubleOptimizedBessel.ComplexBessel.BesselK(ddouble.BitIncrement(nu), x);

                        Console.WriteLine($"{nu}, {x}");
                        Console.WriteLine(y);
                        Console.WriteLine(y_dec);
                        Console.WriteLine(y_inc);

                        Assert.IsTrue(ddouble.Abs((y_dec.R - y.R) / y.R) < 8e-30);
                        Assert.IsTrue(ddouble.Abs((y_inc.R - y.R) / y.R) < 8e-30);
                    }

                    for (ddouble x = 1d / 32; x <= 8; x += 1d / 32) {
                        Complex y = DDoubleOptimizedBessel.ComplexBessel.BesselK(nu, (x, x));
                        Complex y_dec = DDoubleOptimizedBessel.ComplexBessel.BesselK(ddouble.BitDecrement(nu), (x, x));
                        Complex y_inc = DDoubleOptimizedBessel.ComplexBessel.BesselK(ddouble.BitIncrement(nu), (x, x));

                        Console.WriteLine($"{nu}, {(x, x)}");
                        Console.WriteLine(y);
                        Console.WriteLine(y_dec);
                        Console.WriteLine(y_inc);

                        Assert.IsTrue((y_dec - y).Magnitude / y.Magnitude < 8e-30);
                        Assert.IsTrue((y_inc - y).Magnitude / y.Magnitude < 8e-30);
                    }

                    for (ddouble x = 1d / 32; x <= 8; x += 1d / 32) {
                        Complex y = DDoubleOptimizedBessel.ComplexBessel.BesselK(nu, (0, x));
                        Complex y_dec = DDoubleOptimizedBessel.ComplexBessel.BesselK(ddouble.BitDecrement(nu), (0, x));
                        Complex y_inc = DDoubleOptimizedBessel.ComplexBessel.BesselK(ddouble.BitIncrement(nu), (0, x));

                        Console.WriteLine($"{nu}, {(0, x)}");
                        Console.WriteLine(y);
                        Console.WriteLine(y_dec);
                        Console.WriteLine(y_inc);

                        Assert.IsTrue((y_dec - y).Magnitude / y.Magnitude < 8e-30);
                        Assert.IsTrue((y_inc - y).Magnitude / y.Magnitude < 8e-30);
                    }

                    for (ddouble x = 1d / 32; x <= 8; x += 1d / 32) {
                        Complex y = DDoubleOptimizedBessel.ComplexBessel.BesselK(nu, (-x, x));
                        Complex y_dec = DDoubleOptimizedBessel.ComplexBessel.BesselK(ddouble.BitDecrement(nu), (-x, x));
                        Complex y_inc = DDoubleOptimizedBessel.ComplexBessel.BesselK(ddouble.BitIncrement(nu), (-x, x));

                        Console.WriteLine($"{nu}, {(-x, x)}");
                        Console.WriteLine(y);
                        Console.WriteLine(y_dec);
                        Console.WriteLine(y_inc);

                        Assert.IsTrue((y_dec - y).Magnitude / y.Magnitude < 8e-30);
                        Assert.IsTrue((y_inc - y).Magnitude / y.Magnitude < 8e-30);
                    }

                    for (ddouble x = 1d / 4; x <= 8; x += 1d / 4) {
                        Complex y = DDoubleOptimizedBessel.ComplexBessel.BesselK(nu, -x);
                        Complex y_dec = DDoubleOptimizedBessel.ComplexBessel.BesselK(ddouble.BitDecrement(nu), -x);
                        Complex y_inc = DDoubleOptimizedBessel.ComplexBessel.BesselK(ddouble.BitIncrement(nu), -x);

                        Console.WriteLine($"{nu}, {-x}");
                        Console.WriteLine(y);
                        Console.WriteLine(y_dec);
                        Console.WriteLine(y_inc);

                        Assert.IsTrue(ddouble.Abs((y_dec.R - y.R) / y.R) < 8e-30);
                        Assert.IsTrue(ddouble.Abs((y_inc.R - y.R) / y.R) < 8e-30);
                    }
                }
            }

            ddouble eps = Math.ScaleB(1, -9);

            foreach ((ddouble nu, ddouble[] expecteds) in new (ddouble, ddouble[])[] {
                (15 + eps, nu15peps_expecteds), (16 - eps, nu16meps_expecteds),
                (eps, nu0peps_expecteds), (4 - eps, nu4meps_expecteds), (4 + eps, nu4peps_expecteds)
            }) {
                for ((ddouble x, int i) = (1d / 16, 0); i < expecteds.Length; x += 1d / 16, i++) {
                    Complex expected = expecteds[i];

                    Complex actual = DDoubleOptimizedBessel.ComplexBessel.BesselK(nu, x);

                    Console.WriteLine($"{nu}, {x}");
                    Console.WriteLine(expected);
                    Console.WriteLine(actual);

                    Assert.IsTrue(ddouble.Abs((actual.R - expected.R) / expected.R) < 8e-30);
                }
            }
        }
    }
}