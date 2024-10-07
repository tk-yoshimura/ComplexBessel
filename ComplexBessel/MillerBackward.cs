﻿using MultiPrecision;
using MultiPrecisionComplex;
using System.Diagnostics;

namespace ComplexBessel {
    public class MillerBackward<N> where N : struct, IConstant {
        public static readonly double MillerBwdBesselYEps = double.ScaleB(1, -30);

        private static readonly Dictionary<MultiPrecision<N>, BesselJPhiTable> phi_coef_table = [];
        private static readonly Dictionary<MultiPrecision<N>, BesselIPsiTable> psi_coef_table = [];
        private static readonly Dictionary<MultiPrecision<N>, BesselYEtaTable> eta_coef_table = [];
        private static readonly Dictionary<MultiPrecision<N>, BesselYXiTable> xi_coef_table = [];

        public static Complex<N> BesselJ(int n, Complex<N> z) {
            Debug.Assert(MultiPrecision<N>.IsPositive(z.R));
            Debug.Assert(MultiPrecision<N>.IsPositive(z.I));

            int m = BesselJYIterM((double)z.R);

            Complex<N> y = BesselJKernel(n, z, m);

            return y;
        }

        public static Complex<N> BesselJ(MultiPrecision<N> nu, Complex<N> z) {
            Debug.Assert(MultiPrecision<N>.IsPositive(z.R));
            Debug.Assert(MultiPrecision<N>.IsPositive(z.I));

            int m = BesselJYIterM((double)z.R);

            if (BesselUtil<N>.NearlyInteger(nu, out int n)) {
                Complex<N> y = BesselJKernel(n, z, m);

                return y;
            }
            else {
                Complex<N> y = BesselJKernel(nu, z, m);

                return y;
            }
        }

        public static Complex<N> BesselY(int n, Complex<N> z) {
            Debug.Assert(MultiPrecision<N>.IsPositive(z.R));
            Debug.Assert(MultiPrecision<N>.IsPositive(z.I));

            int m = BesselJYIterM((double)z.R);

            Complex<N> y = BesselYKernel(n, z, m);

            return y;
        }

        public static Complex<N> BesselY(MultiPrecision<N> nu, Complex<N> z) {
            Debug.Assert(MultiPrecision<N>.IsPositive(z.R));
            Debug.Assert(MultiPrecision<N>.IsPositive(z.I));

            int m = BesselJYIterM((double)z.R);

            if (BesselUtil<N>.NearlyInteger(nu, out int n)) {
                Complex<N> y = BesselYKernel(n, z, m);

                return y;
            }
            else {
                Complex<N> y = BesselYKernel(nu, z, m);

                return y;
            }
        }

        private static int BesselJYIterM(double r) {
            return 256;
        }

        public static Complex<N> BesselI(int n, Complex<N> z) {
            Debug.Assert(MultiPrecision<N>.IsPositive(z.R));
            Debug.Assert(MultiPrecision<N>.IsPositive(z.I));

            int m = BesselIIterM((double)z.R, (double)z.I);

            Complex<N> y = BesselIKernel(n, z, m);

            return y;
        }

        public static Complex<N> BesselI(MultiPrecision<N> nu, Complex<N> z) {
            Debug.Assert(MultiPrecision<N>.IsPositive(z.R));
            Debug.Assert(MultiPrecision<N>.IsPositive(z.I));

            int m = BesselIIterM((double)z.R, (double)z.I);

            if (BesselUtil<N>.NearlyInteger(nu, out int n)) {
                Complex<N> y = BesselIKernel(n, z, m);

                return y;
            }
            else {
                Complex<N> y = BesselIKernel(nu, z, m);

                return y;
            }
        }

        private static int BesselIIterM(double r, double i) {
            return 256;
        }

        private static Complex<N> BesselJKernel(int n, Complex<N> z, int m) {
            Debug.Assert(m >= 2 && (m & 1) == 0 && n < m);

            if (n < 0) {
                return (n & 1) == 0 ? BesselJKernel(-n, z, m) : -BesselJKernel(-n, z, m);
            }
            else if (n == 0) {
                return BesselJ0Kernel(z, m);
            }
            else if (n == 1) {
                return BesselJ1Kernel(z, m);
            }
            else {
                return BesselJNKernel(n, z, m);
            }
        }

        public static Complex<N> BesselJKernel(MultiPrecision<N> nu, Complex<N> z, int m) {
            int n = (int)MultiPrecision<N>.Floor(nu);
            MultiPrecision<N> alpha = nu - n;

            if (alpha == 0d) {
                return BesselJKernel(n, z, m);
            }

            Debug.Assert(m >= 2 && (m & 1) == 0 && n < m);

            if (!phi_coef_table.TryGetValue(alpha, out BesselJPhiTable phi)) {
                phi = new BesselJPhiTable(alpha);
                phi_coef_table.Add(alpha, phi);
            }

            Complex<N> f0 = 1e-256d, f1 = 0d, lambda = 0d;
            Complex<N> v = 1d / z;

            if (n >= 0) {
                Complex<N> fn = 0d;

                for (int k = m; k >= 1; k--) {
                    if ((k & 1) == 0) {
                        lambda += f0 * phi[k / 2];
                    }

                    (f0, f1) = (Complex<N>.Ldexp(k + alpha, 1) * v * f0 - f1, f0);

                    if (k - 1 == n) {
                        fn = f0;
                    }
                }

                lambda += f0 * phi[0];
                lambda *= Complex<N>.Pow(Complex<N>.Ldexp(v, 1), alpha);

                Complex<N> yn = fn / lambda;

                return yn;
            }
            else {
                for (int k = m; k >= 1; k--) {
                    if ((k & 1) == 0) {
                        lambda += f0 * phi[k / 2];
                    }

                    (f0, f1) = (Complex<N>.Ldexp(k + alpha, 1) * v * f0 - f1, f0);
                }

                lambda += f0 * phi[0];
                lambda *= Complex<N>.Pow(Complex<N>.Ldexp(v, 1), alpha);

                for (int k = 0; k > n; k--) {
                    (f0, f1) = (Complex<N>.Ldexp(k + alpha, 1) * v * f0 - f1, f0);
                }

                Complex<N> yn = f0 / lambda;

                return yn;
            }
        }

        private static Complex<N> BesselJ0Kernel(Complex<N> z, int m) {
            Debug.Assert(m >= 2 && (m & 1) == 0);

            Complex<N> f0 = 1e-256d, f1 = 0d, lambda = 0d;
            Complex<N> v = 1d / z;

            for (int k = m; k >= 1; k--) {
                if ((k & 1) == 0) {
                    lambda += f0;
                }

                (f0, f1) = (2 * k * v * f0 - f1, f0);
            }

            lambda = Complex<N>.Ldexp(lambda, 1) + f0;

            Complex<N> y0 = f0 / lambda;

            return y0;
        }

        private static Complex<N> BesselJ1Kernel(Complex<N> z, int m) {
            Debug.Assert(m >= 2 && (m & 1) == 0);

            Complex<N> f0 = 1e-256d, f1 = 0d, lambda = 0d;
            Complex<N> v = 1d / z;

            for (int k = m; k >= 1; k--) {
                if ((k & 1) == 0) {
                    lambda += f0;
                }

                (f0, f1) = (2 * k * v * f0 - f1, f0);
            }

            lambda = Complex<N>.Ldexp(lambda, 1) + f0;

            Complex<N> y1 = f1 / lambda;

            return y1;
        }

        private static Complex<N> BesselJNKernel(int n, Complex<N> z, int m) {
            Complex<N> f0 = 1e-256d, f1 = 0d, fn = 0d, lambda = 0d;
            Complex<N> v = 1d / z;

            for (int k = m; k >= 1; k--) {
                if ((k & 1) == 0) {
                    lambda += f0;
                }

                (f0, f1) = (2 * k * v * f0 - f1, f0);

                if (k - 1 == n) {
                    fn = f0;
                }
            }

            lambda = Complex<N>.Ldexp(lambda, 1) + f0;

            Complex<N> yn = fn / lambda;

            return yn;
        }

        private static Complex<N> BesselYKernel(int n, Complex<N> z, int m) {
            Debug.Assert(m >= 2 && (m & 1) == 0 && n < m);

            if (n < 0) {
                return (n & 1) == 0 ? BesselYKernel(-n, z, m) : -BesselYKernel(-n, z, m);
            }
            else if (n == 0) {
                return BesselY0Kernel(z, m);
            }
            else if (n == 1) {
                return BesselY1Kernel(z, m);
            }
            else {
                return BesselYNKernel(n, z, m);
            }
        }

        public static Complex<N> BesselYKernel(MultiPrecision<N> nu, Complex<N> z, int m) {
            int n = (int)MultiPrecision<N>.Floor(nu);
            MultiPrecision<N> alpha = nu - n;

            if (alpha == 0d) {
                return BesselYKernel(n, z, m);
            }

            Debug.Assert(m >= 2 && (m & 1) == 0 && n < m);

            if (!eta_coef_table.TryGetValue(alpha, out BesselYEtaTable eta)) {
                eta = new BesselYEtaTable(alpha);
                eta_coef_table.Add(alpha, eta);
            }
            if (!xi_coef_table.TryGetValue(alpha, out BesselYXiTable xi)) {
                xi = new BesselYXiTable(alpha, eta);
                xi_coef_table.Add(alpha, xi);
            }
            if (!phi_coef_table.TryGetValue(alpha, out BesselJPhiTable phi)) {
                phi = new BesselJPhiTable(alpha);
                phi_coef_table.Add(alpha, phi);
            }

            Complex<N> f0 = 1e-256, f1 = 0d, lambda = 0d;
            Complex<N> se = 0d, sxo = 0d, sxe = 0d;
            Complex<N> v = 1d / z;

            for (int k = m; k >= 1; k--) {
                if ((k & 1) == 0) {
                    lambda += f0 * phi[k / 2];

                    se += f0 * eta[k / 2];
                    sxe += f0 * xi[k];
                }
                else if (k >= 3) {
                    sxo += f0 * xi[k];
                }

                (f0, f1) = (Complex<N>.Ldexp(k + alpha, 1) * v * f0 - f1, f0);
            }

            Complex<N> s = Complex<N>.Pow(Complex<N>.Ldexp(v, 1), alpha), sqs = s * s;

            lambda += f0 * phi[0];
            lambda *= s;

            MultiPrecision<N> rcot = 1d / MultiPrecision<N>.TanPI(alpha), rgamma = MultiPrecision<N>.Gamma(1d + alpha), rsqgamma = rgamma * rgamma;
            Complex<N> r = Complex<N>.Ldexp(MultiPrecision<N>.RcpPI * sqs, 1);
            Complex<N> p = sqs * rsqgamma * MultiPrecision<N>.RcpPI;

            Complex<N> eta0 = MultiPrecision<N>.Abs(alpha) > MillerBwdBesselYEps
                ? rcot - p / alpha
                : BesselYEta0Eps(alpha, z);

            Complex<N> xi0 = -Complex<N>.Ldexp(v, 1) * p;
            Complex<N> xi1 = MultiPrecision<N>.Abs(alpha) > MillerBwdBesselYEps
                ? rcot + p * (alpha * (alpha + 1d) + 1d) / (alpha * (alpha - 1d))
                : BesselYXi1Eps(alpha, z);

            Complex<N> y0 = r * se + eta0 * f0;
            Complex<N> y1 = r * (3d * alpha * v * sxe + sxo) + xi0 * f0 + xi1 * f1;

            if (n == 0) {
                Complex<N> yn = y0 / lambda;

                return yn;
            }
            if (n == 1) {
                Complex<N> yn = y1 / lambda;

                return yn;
            }
            if (n >= 0) {
                for (int k = 1; k < n; k++) {
                    (y1, y0) = (Complex<N>.Ldexp(k + alpha, 1) * v * y1 - y0, y1);
                }

                Complex<N> yn = y1 / lambda;

                return yn;
            }
            else {
                for (int k = 0; k > n; k--) {
                    (y0, y1) = (Complex<N>.Ldexp(k + alpha, 1) * v * y0 - y1, y0);
                }

                Complex<N> yn = y0 / lambda;

                return yn;
            }
        }

        private static Complex<N> BesselYEta0Eps(MultiPrecision<N> alpha, Complex<N> z) {
            Complex<N> lnz = Complex<N>.Log(z), lnhalfz = Complex<N>.Log(Complex<N>.Ldexp(z, -1));
            MultiPrecision<N> pi = MultiPrecision<N>.PI, sqpi = pi * pi;
            MultiPrecision<N> ln2 = MultiPrecision<N>.Ln2, sqln2 = ln2 * ln2, cbln2 = sqln2 * ln2, qdln2 = sqln2 * sqln2;
            MultiPrecision<N> g = MultiPrecision<N>.EulerGamma;

            Complex<N> r0 = lnhalfz + g;
            Complex<N> r1 =
                (-sqln2 + lnz * (ln2 * 2d - lnz)) * 4d
                - sqpi
                - g * (lnhalfz * 2d + g) * 4d;
            Complex<N> r2 =
                (-cbln2 + lnz * (sqln2 * 3d + lnz * (ln2 * -3d + lnz))) * 4d
                + MultiPrecision<N>.Zeta3 * 2d
                + sqpi * (lnhalfz + g)
                + g * ((sqln2 + lnz * (ln2 * -2d + lnz)) * 3d + g * (lnhalfz * 3d + g)) * 4d;
            Complex<N> r3 =
                (-qdln2 + lnz * (cbln2 * 4d + lnz * (sqln2 * -6d + lnz * (ln2 * 4d - lnz)))) * 16d
                - MultiPrecision<N>.Zeta3 * (lnhalfz + g) * 32d
                - sqpi * ((sqln2 + lnz * (-ln2 * 2d + lnz) + g * (lnhalfz * 2d + g)) * 8d + sqpi)
                + g * ((cbln2 + lnz * (sqln2 * -3d + lnz * (ln2 * 3d - lnz))) * 4d
                + g * ((sqln2 + lnz * (ln2 * -2d + lnz)) * -6d
                + g * (lnhalfz * -4d
                - g))) * 16d;

            Complex<N> eta0 = (r0 * 48d + alpha * (r1 * 12d + alpha * (r2 * 8d + alpha * r3))) / (24d * MultiPrecision<N>.PI);

            return eta0;
        }

        private static Complex<N> BesselYXi1Eps(MultiPrecision<N> alpha, Complex<N> z) {
            Complex<N> lnz = Complex<N>.Log(z), lnhalfz = Complex<N>.Log(Complex<N>.Ldexp(z, -1)), lnzm1 = lnz - 1, lnhalfzm1 = lnhalfz - 1;
            MultiPrecision<N> pi = MultiPrecision<N>.PI, sqpi = pi * pi;
            MultiPrecision<N> ln2 = MultiPrecision<N>.Ln2, sqln2 = ln2 * ln2, cbln2 = sqln2 * ln2, qdln2 = sqln2 * sqln2;
            MultiPrecision<N> g = MultiPrecision<N>.EulerGamma;

            Complex<N> r0 = lnhalfzm1 + g;
            Complex<N> r1 =
                (-sqln2 + ln2 * lnzm1 * 2d + lnz * (2 - lnz)) * 4d
                - sqpi
                - g * (lnhalfzm1 * 2d + g) * 4d
                - 6d;
            Complex<N> r2 =
                -cbln2 * 4d + sqln2 * lnzm1 * 12d + lnz * (18d + lnz * (-12d + lnz * 4d))
                + ln2 * (lnz * (2d - lnz) * 12d - 18d)
                + MultiPrecision<N>.Zeta3 * 2d
                + sqpi * (lnhalfzm1 + g)
                + g * ((sqln2 - ln2 * lnzm1 * 2d + lnz * (-2d + lnz)) * 12d + 18d
                + g * (lnhalfzm1 * 12d
                + g * 4d))
                - 9d;
            Complex<N> r3 =
                -qdln2 * 16d
                + cbln2 * lnzm1 * 64d
                + sqln2 * (lnz * (2d - lnz) * 96d - 144d)
                + ln2 * (lnz * (9d + lnz * (-6d + lnz * 2d)) * 32d - 144d)
                + lnz * (9d + lnz * (-9d + lnz * (4d - lnz))) * 16d
                + MultiPrecision<N>.Zeta3 * (lnhalfzm1 + g) * -32d
                + sqpi * ((-sqln2 + ln2 * lnzm1 * 2d + lnz * (2d - lnz) - g * (lnhalfzm1 * 2d + g)) * 8d - 12d - sqpi)
                + g * ((cbln2 - sqln2 * lnzm1 * 3d) * 64d + ln2 * (lnz * (-2d + lnz) * 192d + 288d) + lnz * (-9d + lnz * (6d - lnz * 2d)) * 32d + 144d
                + g * ((-sqln2 + ln2 * lnzm1 * 2d + lnz * (2d - lnz)) * 96d - 144d
                + g * (lnhalfzm1 * -64d
                - g * 16d)))
                - 72d;

            Complex<N> xi1 = (r0 * 48d + alpha * (r1 * 12d + alpha * (r2 * 8d + alpha * r3))) / (24d * MultiPrecision<N>.PI);

            return xi1;
        }

        private static Complex<N> BesselY0Kernel(Complex<N> z, int m) {
            Debug.Assert(m >= 2 && (m & 1) == 0);

            if (!eta_coef_table.TryGetValue(0, out BesselYEtaTable eta)) {
                eta = new BesselYEtaTable(0);
                eta_coef_table.Add(0, eta);
            }

            Complex<N> f0 = 1e-256, f1 = 0d, lambda = 0d;
            Complex<N> se = 0d;
            Complex<N> v = 1d / z;

            for (int k = m; k >= 1; k--) {
                if ((k & 1) == 0) {
                    lambda += f0;

                    se += f0 * eta[k / 2];
                }

                (f0, f1) = (2 * k * v * f0 - f1, f0);
            }

            lambda = Complex<N>.Ldexp(lambda, 1) + f0;

            Complex<N> y0 = Complex<N>.Ldexp((se + f0 * (Complex<N>.Log(Complex<N>.Ldexp(z, -1)) + MultiPrecision<N>.EulerGamma)) / (MultiPrecision<N>.PI * lambda), 1);

            return y0;
        }

        private static Complex<N> BesselY1Kernel(Complex<N> z, int m) {
            Debug.Assert(m >= 2 && (m & 1) == 0);

            if (!xi_coef_table.TryGetValue(0, out BesselYXiTable xi)) {
                if (!eta_coef_table.TryGetValue(0, out BesselYEtaTable eta)) {
                    eta = new BesselYEtaTable(0);
                    eta_coef_table.Add(0, eta);
                }

                xi = new BesselYXiTable(0, eta);
                xi_coef_table.Add(0, xi);
            }

            Complex<N> f0 = 1e-256, f1 = 0d, lambda = 0d;
            Complex<N> sx = 0d;
            Complex<N> v = 1d / z;

            for (int k = m; k >= 1; k--) {
                if ((k & 1) == 0) {
                    lambda += f0;
                }
                else if (k >= 3) {
                    sx += f0 * xi[k];
                }

                (f0, f1) = (2 * k * v * f0 - f1, f0);
            }

            lambda = Complex<N>.Ldexp(lambda, 1) + f0;

            Complex<N> y1 = Complex<N>.Ldexp((sx - v * f0 + (Complex<N>.Log(Complex<N>.Ldexp(z, -1)) + MultiPrecision<N>.EulerGamma - 1d) * f1) / (lambda * MultiPrecision<N>.PI), 1);

            return y1;
        }

        private static Complex<N> BesselYNKernel(int n, Complex<N> z, int m) {
            if (!eta_coef_table.TryGetValue(0, out BesselYEtaTable eta)) {
                eta = new BesselYEtaTable(0);
                eta_coef_table.Add(0, eta);
            }

            if (!xi_coef_table.TryGetValue(0, out BesselYXiTable xi)) {
                xi = new BesselYXiTable(0, eta);
                xi_coef_table.Add(0, xi);
            }

            Complex<N> f0 = 1e-256, f1 = 0d, lambda = 0d;
            Complex<N> se = 0d, sx = 0d;
            Complex<N> v = 1d / z;

            for (int k = m; k >= 1; k--) {
                if ((k & 1) == 0) {
                    lambda += f0;

                    se += f0 * eta[k / 2];
                }
                else if (k >= 3) {
                    sx += f0 * xi[k];
                }

                (f0, f1) = (2 * k * v * f0 - f1, f0);
            }

            lambda = Complex<N>.Ldexp(lambda, 1) + f0;

            Complex<N> c = Complex<N>.Log(Complex<N>.Ldexp(z, -1)) + MultiPrecision<N>.EulerGamma;

            Complex<N> y0 = se + f0 * c;
            Complex<N> y1 = sx - v * f0 + (c - 1d) * f1;

            for (int k = 1; k < n; k++) {
                (y1, y0) = (2 * k * v * y1 - y0, y1);
            }

            Complex<N> yn = Complex<N>.Ldexp(y1 / (lambda * MultiPrecision<N>.PI), 1);

            return yn;
        }

        private static Complex<N> BesselIKernel(int n, Complex<N> z, int m) {
            Debug.Assert(m >= 2 && (m & 1) == 0 && n < m);

            n = int.Abs(n);

            if (n == 0) {
                return BesselI0Kernel(z, m);
            }
            else if (n == 1) {
                return BesselI1Kernel(z, m);
            }
            else {
                return BesselINKernel(n, z, m);
            }
        }

        public static Complex<N> BesselIKernel(MultiPrecision<N> nu, Complex<N> z, int m) {
            int n = (int)MultiPrecision<N>.Floor(nu);
            MultiPrecision<N> alpha = nu - n;

            if (alpha == 0d) {
                return BesselIKernel(n, z, m);
            }

            Debug.Assert(m >= 2 && (m & 1) == 0 && n < m);

            if (!psi_coef_table.TryGetValue(alpha, out BesselIPsiTable psi)) {
                psi = new BesselIPsiTable(alpha);
                psi_coef_table.Add(alpha, psi);
            }

            Complex<N> g0 = 1e-256, g1 = 0d, lambda = 0d;
            Complex<N> v = 1d / z;

            if (n >= 0) {
                Complex<N> gn = 0d;

                for (int k = m; k >= 1; k--) {
                    lambda += g0 * psi[k];

                    (g0, g1) = (Complex<N>.Ldexp(k + alpha, 1) * v * g0 + g1, g0);

                    if (k - 1 == n) {
                        gn = g0;
                    }
                }

                lambda += g0 * psi[0];
                lambda *= Complex<N>.Pow(Complex<N>.Ldexp(v, 1), alpha);

                Complex<N> yn = gn / lambda;

                yn *= Complex<N>.Exp(z);

                return yn;
            }
            else {
                for (int k = m; k >= 1; k--) {
                    lambda += g0 * psi[k];

                    (g0, g1) = (Complex<N>.Ldexp(k + alpha, 1) * v * g0 + g1, g0);
                }

                lambda += g0 * psi[0];
                lambda *= Complex<N>.Pow(Complex<N>.Ldexp(v, 1), alpha);

                for (int k = 0; k > n; k--) {
                    (g0, g1) = (Complex<N>.Ldexp(k + alpha, 1) * v * g0 + g1, g0);
                }

                Complex<N> yn = g0 / lambda;

                yn *= Complex<N>.Exp(z);

                return yn;
            }
        }

        private static Complex<N> BesselI0Kernel(Complex<N> z, int m) {
            Debug.Assert(m >= 2 && (m & 1) == 0);

            Complex<N> g0 = 1e-256, g1 = 0d, lambda = 0d;
            Complex<N> v = 1d / z;

            for (int k = m; k >= 1; k--) {
                lambda += g0;

                (g0, g1) = (2 * k * v * g0 + g1, g0);
            }

            lambda = Complex<N>.Ldexp(lambda, 1) + g0;

            Complex<N> y0 = g0 / lambda;

            y0 *= Complex<N>.Exp(z);

            return y0;
        }

        private static Complex<N> BesselI1Kernel(Complex<N> z, int m) {
            Debug.Assert(m >= 2 && (m & 1) == 0);

            Complex<N> g0 = 1e-256, g1 = 0d, lambda = 0d;
            Complex<N> v = 1d / z;

            for (int k = m; k >= 1; k--) {
                lambda += g0;

                (g0, g1) = (2 * k * v * g0 + g1, g0);
            }

            lambda = Complex<N>.Ldexp(lambda, 1) + g0;

            Complex<N> y1 = g1 / lambda;

            y1 *= Complex<N>.Exp(z);

            return y1;
        }

        private static Complex<N> BesselINKernel(int n, Complex<N> z, int m) {
            Complex<N> f0 = 1e-256, f1 = 0d, lambda = 0d, fn = 0d;
            Complex<N> v = 1d / z;

            for (int k = m; k >= 1; k--) {
                lambda += f0;

                (f0, f1) = (2 * k * v * f0 + f1, f0);

                if (k - 1 == n) {
                    fn = f0;
                }
            }

            lambda = Complex<N>.Ldexp(lambda, 1) + f0;

            Complex<N> yn = fn / lambda;

            yn *= Complex<N>.Exp(z);

            return yn;
        }

        private class BesselJPhiTable {
            private readonly MultiPrecision<N> alpha;
            private readonly List<MultiPrecision<N>> table = [];

            private MultiPrecision<N> g;

            public BesselJPhiTable(MultiPrecision<N> alpha) {
                Debug.Assert(alpha > 0d && alpha < 1d, nameof(alpha));

                this.alpha = alpha;

                MultiPrecision<N> phi0 = MultiPrecision<N>.Gamma(1 + alpha);
                MultiPrecision<N> phi1 = phi0 * (alpha + 2d);

                this.g = phi0;

                this.table.Add(phi0);
                this.table.Add(phi1);
            }

            public MultiPrecision<N> this[int k] => Value(k);

            private MultiPrecision<N> Value(int k) {
                Debug.Assert(k >= 0);

                if (k < table.Count) {
                    return table[k];
                }

                for (int i = table.Count; i <= k; i++) {
                    g = g * (alpha + i - 1d) / i;

                    MultiPrecision<N> phi = g * (alpha + 2 * i);

                    table.Add(phi);
                }

                return table[k];
            }
        };

        private class BesselIPsiTable {
            private readonly MultiPrecision<N> alpha;
            private readonly List<MultiPrecision<N>> table = [];

            private MultiPrecision<N> g;

            public BesselIPsiTable(MultiPrecision<N> alpha) {
                Debug.Assert(alpha > 0d && alpha < 1d, nameof(alpha));

                this.alpha = alpha;

                MultiPrecision<N> psi0 = MultiPrecision<N>.Gamma(1d + alpha);
                MultiPrecision<N> psi1 = MultiPrecision<N>.Ldexp(psi0, 1) * (1d + alpha);

                this.g = MultiPrecision<N>.Ldexp(psi0, 1);

                this.table.Add(psi0);
                this.table.Add(psi1);
            }

            public MultiPrecision<N> this[int k] => Value(k);

            private MultiPrecision<N> Value(int k) {
                Debug.Assert(k >= 0);

                if (k < table.Count) {
                    return table[k];
                }

                for (int i = table.Count; i <= k; i++) {
                    g = g * (MultiPrecision<N>.Ldexp(alpha, 1) + i - 1d) / i;

                    MultiPrecision<N> phi = g * (alpha + i);

                    table.Add(phi);
                }

                return table[k];
            }
        };

        private class BesselYEtaTable {
            private readonly MultiPrecision<N> alpha;
            private readonly List<MultiPrecision<N>> table = [];

            private MultiPrecision<N> g;

            public BesselYEtaTable(MultiPrecision<N> alpha) {
                Debug.Assert(alpha >= 0d && alpha < 1d, nameof(alpha));

                this.alpha = alpha;
                this.table.Add(MultiPrecision<N>.NaN);

                if (alpha > 0d) {
                    MultiPrecision<N> c = MultiPrecision<N>.Gamma(1d + alpha);
                    c *= c;
                    this.g = 1d / (1d - alpha) * c;

                    MultiPrecision<N> eta1 = (alpha + 2d) * g;

                    this.table.Add(eta1);
                }
            }

            public MultiPrecision<N> this[int k] => Value(k);

            private MultiPrecision<N> Value(int k) {
                Debug.Assert(k >= 0);

                if (k < table.Count) {
                    return table[k];
                }

                for (int i = table.Count; i <= k; i++) {
                    if (alpha > 0d) {
                        g = -g * (alpha + i - 1) * (MultiPrecision<N>.Ldexp(alpha, 1) + i - 1d) / (i * (i - alpha));

                        MultiPrecision<N> eta = g * (alpha + 2 * i);

                        table.Add(eta);
                    }
                    else {
                        MultiPrecision<N> eta = (MultiPrecision<N>)2d / i;

                        table.Add((i & 1) == 1 ? eta : -eta);
                    }
                }

                return table[k];
            }
        };

        private class BesselYXiTable {
            private readonly MultiPrecision<N> alpha;
            private readonly List<MultiPrecision<N>> table = [];
            private readonly BesselYEtaTable eta;

            public BesselYXiTable(MultiPrecision<N> alpha, BesselYEtaTable eta) {
                Debug.Assert(alpha >= 0d && alpha < 1d, nameof(alpha));

                this.alpha = alpha;
                this.table.Add(MultiPrecision<N>.NaN);
                this.table.Add(MultiPrecision<N>.NaN);

                this.eta = eta;
            }

            public MultiPrecision<N> this[int k] => Value(k);

            private MultiPrecision<N> Value(int k) {
                Debug.Assert(k >= 0);

                if (k < table.Count) {
                    return table[k];
                }

                for (int i = table.Count; i <= k; i++) {
                    if (alpha > 0d) {
                        if ((i & 1) == 0) {
                            table.Add(eta[i / 2]);
                        }
                        else {
                            table.Add((eta[i / 2] - eta[i / 2 + 1]) / 2);
                        }
                    }
                    else {
                        if ((i & 1) == 1) {
                            MultiPrecision<N> xi = (MultiPrecision<N>)(2 * (i / 2) + 1) / (i / 2 * (i / 2 + 1));
                            table.Add((i & 2) > 0 ? xi : -xi);
                        }
                        else {
                            table.Add(MultiPrecision<N>.NaN);
                        }
                    }
                }

                return table[k];
            }
        }
    }
}
