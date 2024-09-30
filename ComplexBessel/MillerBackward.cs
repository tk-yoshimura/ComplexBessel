using MultiPrecision;
using MultiPrecisionComplex;
using System.Diagnostics;

namespace ComplexBessel {
    public class MillerBackward<N> where N : struct, IConstant {
        private static readonly Dictionary<MultiPrecision<N>, BesselJPhiTable> phi_coef_table = [];
        private static readonly Dictionary<MultiPrecision<N>, BesselIPsiTable> psi_coef_table = [];
        private static readonly Dictionary<MultiPrecision<N>, BesselYEtaTable> eta_coef_table = [];
        private static readonly Dictionary<MultiPrecision<N>, BesselYXiTable> xi_coef_table = [];

        public static Complex<N> BesselJ(int n, Complex<N> x) {
            int m = BesselJYIterM((double)x.Magnitude);

            Complex<N> y = BesselJKernel(n, x, m);

            return y;
        }

        public static Complex<N> BesselJ(MultiPrecision<N> nu, Complex<N> x) {
            int m = BesselJYIterM((double)x.Magnitude);

            if (BesselUtil<N>.NearlyInteger(nu, out int n)) {
                Complex<N> y = BesselJKernel(n, x, m);

                return y;
            }
            else {
                Complex<N> y = BesselJKernel(nu, x, m);

                return y;
            }
        }

        public static Complex<N> BesselY(int n, Complex<N> x) {
            int m = BesselJYIterM((double)x.Magnitude);

            Complex<N> y = BesselYKernel(n, x, m);

            return y;
        }

        public static Complex<N> BesselY(MultiPrecision<N> nu, Complex<N> x) {
            int m = BesselJYIterM((double)x.Magnitude);

            if (BesselUtil<N>.NearlyInteger(nu, out int n)) {
                Complex<N> y = BesselYKernel(n, x, m);

                return y;
            }
            else {
                Complex<N> y = BesselYKernel(nu, x, m);

                return y;
            }
        }

        private static int BesselJYIterM(double x) {
            return 256;//(int)double.Ceiling(74 + 1.36 * x - 54.25 / double.Sqrt(double.Sqrt(x))) & ~1;
        }

        public static Complex<N> BesselI(int n, Complex<N> x, bool scale = false) {
            int m = BesselIIterM((double)x.Magnitude);

            Complex<N> y = BesselIKernel(n, x, m, scale);

            return y;
        }

        public static Complex<N> BesselI(MultiPrecision<N> nu, Complex<N> x, bool scale = false) {
            int m = BesselIIterM((double)x.Magnitude);

            if (BesselUtil<N>.NearlyInteger(nu, out int n)) {
                Complex<N> y = BesselIKernel(n, x, m, scale);

                return y;
            }
            else {
                Complex<N> y = BesselIKernel(nu, x, m, scale);

                return y;
            }
        }

        private static int BesselIIterM(double x) {
            return 256; // (int)double.Ceiling(86 + 0.75 * x - 67.25 / double.Sqrt(double.Sqrt(x))) & ~1;
        }

        internal static Complex<N> BesselJKernel(int n, Complex<N> x, int m) {
            if (m < 2 || (m & 1) != 0 || n >= m) {
                throw new ArgumentOutOfRangeException(nameof(m));
            }

            if (n < 0) {
                return ((n & 1) == 0) ? BesselJKernel(-n, x, m) : -BesselJKernel(-n, x, m);
            }
            if (n == 0) {
                return BesselJ0Kernel(x, m);
            }
            if (n == 1) {
                return BesselJ1Kernel(x, m);
            }

            Complex<N> f0 = 1e-256d, f1 = 0d, fn = 0d, lambda = 0d;
            Complex<N> v = 1d / x;

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

        internal static Complex<N> BesselJKernel(MultiPrecision<N> nu, Complex<N> x, int m) {
            int n = (int)MultiPrecision<N>.Floor(nu);
            MultiPrecision<N> alpha = nu - n;

            if (alpha == 0d) {
                return BesselJKernel(n, x, m);
            }

            if (m < 2 || (m & 1) != 0 || n >= m) {
                throw new ArgumentOutOfRangeException(nameof(m));
            }

            if (!phi_coef_table.TryGetValue(alpha, out BesselJPhiTable phi_table)) {
                phi_table = new BesselJPhiTable(alpha);
                phi_coef_table.Add(alpha, phi_table);
            }

            BesselJPhiTable phi = phi_table;

            Complex<N> f0 = 1e-256d, f1 = 0d, lambda = 0d;
            Complex<N> v = 1d / x;

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

        private static Complex<N> BesselJ0Kernel(Complex<N> x, int m) {
            if (m < 2 || (m & 1) != 0) {
                throw new ArgumentOutOfRangeException(nameof(m));
            }

            Complex<N> f0 = 1e-256d, f1 = 0d, lambda = 0d;
            Complex<N> v = 1d / x;

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

        private static Complex<N> BesselJ1Kernel(Complex<N> x, int m) {
            if (m < 2 || (m & 1) != 0) {
                throw new ArgumentOutOfRangeException(nameof(m));
            }

            Complex<N> f0 = 1e-256d, f1 = 0d, lambda = 0d;
            Complex<N> v = 1d / x;

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

        internal static Complex<N> BesselYKernel(int n, Complex<N> x, int m) {
            if (m < 2 || (m & 1) != 0 || n >= m) {
                throw new ArgumentOutOfRangeException(nameof(m));
            }

            if (n < 0) {
                return ((n & 1) == 0) ? BesselYKernel(-n, x, m) : -BesselYKernel(-n, x, m);
            }
            if (n == 0) {
                return BesselY0Kernel(x, m);
            }
            if (n == 1) {
                return BesselY1Kernel(x, m);
            }

            if (!eta_coef_table.ContainsKey(0)) {
                eta_coef_table.Add(0, new BesselYEtaTable(0));
            }

            BesselYEtaTable eta = eta_coef_table[0];

            if (!xi_coef_table.ContainsKey(0)) {
                xi_coef_table.Add(0, new BesselYXiTable(0, eta));
            }

            BesselYXiTable xi = xi_coef_table[0];

            Complex<N> f0 = 1e-256, f1 = 0d, lambda = 0d;
            Complex<N> se = 0d, sx = 0d;
            Complex<N> v = 1d / x;

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

            Complex<N> c = Complex<N>.Log(Complex<N>.Ldexp(x, -1)) + MultiPrecision<N>.EulerGamma;

            Complex<N> y0 = se + f0 * c;
            Complex<N> y1 = sx - v * f0 + (c - 1d) * f1;

            for (int k = 1; k < n; k++) {
                (y1, y0) = (2 * k * v * y1 - y0, y1);
            }

            Complex<N> yn = Complex<N>.Ldexp(y1 / (lambda * MultiPrecision<N>.PI), 1);

            return yn;
        }

        internal static Complex<N> BesselYKernel(MultiPrecision<N> nu, Complex<N> x, int m) {
            int n = (int)MultiPrecision<N>.Floor(nu);
            MultiPrecision<N> alpha = nu - n;

            if (alpha == 0d) {
                return BesselYKernel(n, x, m);
            }

            if (m < 2 || (m & 1) != 0 || n >= m) {
                throw new ArgumentOutOfRangeException(nameof(m));
            }

            if (!eta_coef_table.TryGetValue(alpha, out BesselYEtaTable eta_table)) {
                eta_table = new BesselYEtaTable(alpha);
                eta_coef_table.Add(alpha, eta_table);
            }

            BesselYEtaTable eta = eta_table;

            if (!xi_coef_table.TryGetValue(alpha, out BesselYXiTable xi_table)) {
                xi_table = new BesselYXiTable(alpha, eta);
                xi_coef_table.Add(alpha, xi_table);
            }

            BesselYXiTable xi = xi_table;

            if (!phi_coef_table.TryGetValue(alpha, out BesselJPhiTable phi_table)) {
                phi_table = new BesselJPhiTable(alpha);
                phi_coef_table.Add(alpha, phi_table);
            }

            BesselJPhiTable phi = phi_table;

            Complex<N> f0 = 1e-256, f1 = 0d, lambda = 0d;
            Complex<N> se = 0d, sxo = 0d, sxe = 0d;
            Complex<N> v = 1d / x;

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

            Complex<N> eta0 = (MultiPrecision<N>.Abs(alpha) > BesselUtil<N>.MillerBwdBesselYEps)
                ? (rcot - p / alpha)
                : BesselYEta0Eps(alpha, x);

            Complex<N> xi0 = -Complex<N>.Ldexp(v, 1) * p;
            Complex<N> xi1 = (MultiPrecision<N>.Abs(alpha) > BesselUtil<N>.MillerBwdBesselYEps)
                ? rcot + p * (alpha * (alpha + 1d) + 1d) / (alpha * (alpha - 1d))
                : BesselYXi1Eps(alpha, x);

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

        private static Complex<N> BesselYEta0Eps(MultiPrecision<N> alpha, Complex<N> x) {
            Complex<N> lnx = Complex<N>.Log(x), lnhalfx = Complex<N>.Log(Complex<N>.Ldexp(x, -1));
            MultiPrecision<N> pi = MultiPrecision<N>.PI, sqpi = pi * pi;
            MultiPrecision<N> ln2 = MultiPrecision<N>.Ln2, sqln2 = ln2 * ln2, cbln2 = sqln2 * ln2, qdln2 = sqln2 * sqln2;
            MultiPrecision<N> g = MultiPrecision<N>.EulerGamma;

            Complex<N> r0 = lnhalfx + g;
            Complex<N> r1 =
                (-sqln2 + lnx * (ln2 * 2d - lnx)) * 4d
                - sqpi
                - g * (lnhalfx * 2d + g) * 4d;
            Complex<N> r2 =
                (-cbln2 + lnx * (sqln2 * 3d + lnx * (ln2 * -3d + lnx))) * 4d
                + MultiPrecision<N>.Zeta3 * 2d
                + sqpi * (lnhalfx + g)
                + g * ((sqln2 + lnx * (ln2 * -2d + lnx)) * 3d + g * (lnhalfx * 3d + g)) * 4d;
            Complex<N> r3 =
                (-qdln2 + lnx * (cbln2 * 4d + lnx * (sqln2 * -6d + lnx * (ln2 * 4d - lnx)))) * 16d
                - MultiPrecision<N>.Zeta3 * (lnhalfx + g) * 32d
                - sqpi * (((sqln2 + lnx * (-ln2 * 2d + lnx) + g * (lnhalfx * 2d + g)) * 8d) + sqpi)
                + g * ((cbln2 + lnx * (sqln2 * -3d + lnx * (ln2 * 3d - lnx))) * 4d
                + g * ((sqln2 + lnx * (ln2 * -2d + lnx)) * -6d
                + g * (lnhalfx * -4d
                - g))) * 16d;

            Complex<N> eta0 = (r0 * 48d + alpha * (r1 * 12d + alpha * (r2 * 8d + alpha * r3))) / (24d * MultiPrecision<N>.PI);

            return eta0;
        }

        static Complex<N> BesselYXi1Eps(MultiPrecision<N> alpha, Complex<N> x) {
            Complex<N> lnx = Complex<N>.Log(x), lnhalfx = Complex<N>.Log(Complex<N>.Ldexp(x, -1)), lnxm1 = lnx - 1, lnhalfxm1 = lnhalfx - 1;
            MultiPrecision<N> pi = MultiPrecision<N>.PI, sqpi = pi * pi;
            MultiPrecision<N> ln2 = MultiPrecision<N>.Ln2, sqln2 = ln2 * ln2, cbln2 = sqln2 * ln2, qdln2 = sqln2 * sqln2;
            MultiPrecision<N> g = MultiPrecision<N>.EulerGamma;

            Complex<N> r0 = lnhalfxm1 + g;
            Complex<N> r1 =
                (-sqln2 + ln2 * lnxm1 * 2d + lnx * (2 - lnx)) * 4d
                - sqpi
                - g * (lnhalfxm1 * 2d + g) * 4d
                - 6d;
            Complex<N> r2 =
                -cbln2 * 4d + sqln2 * lnxm1 * 12d + lnx * (18d + lnx * (-12d + lnx * 4d))
                + ln2 * (lnx * (2d - lnx) * 12d - 18d)
                + MultiPrecision<N>.Zeta3 * 2d
                + sqpi * (lnhalfxm1 + g)
                + g * ((sqln2 - ln2 * lnxm1 * 2d + lnx * (-2d + lnx)) * 12d + 18d
                + g * (lnhalfxm1 * 12d
                + g * 4d))
                - 9d;
            Complex<N> r3 =
                -qdln2 * 16d
                + cbln2 * lnxm1 * 64d
                + sqln2 * (lnx * (2d - lnx) * 96d - 144d)
                + ln2 * (lnx * (9d + lnx * (-6d + lnx * 2d)) * 32d - 144d)
                + lnx * (9d + lnx * (-9d + lnx * (4d - lnx))) * 16d
                + MultiPrecision<N>.Zeta3 * (lnhalfxm1 + g) * -32d
                + sqpi * ((-sqln2 + ln2 * lnxm1 * 2d + lnx * (2d - lnx) - g * (lnhalfxm1 * 2d + g)) * 8d - 12d - sqpi)
                + g * ((cbln2 - sqln2 * lnxm1 * 3d) * 64d + ln2 * (lnx * (-2d + lnx) * 192d + 288d) + lnx * (-9d + lnx * (6d - lnx * 2d)) * 32d + 144d
                + g * ((-sqln2 + ln2 * lnxm1 * 2d + lnx * (2d - lnx)) * 96d - 144d
                + g * (lnhalfxm1 * -64d
                - g * 16d)))
                - 72d;

            Complex<N> xi1 = (r0 * 48d + alpha * (r1 * 12d + alpha * (r2 * 8d + alpha * r3))) / (24d * MultiPrecision<N>.PI);

            return xi1;
        }

        private static Complex<N> BesselY0Kernel(Complex<N> x, int m) {
            if (m < 2 || (m & 1) != 0) {
                throw new ArgumentOutOfRangeException(nameof(m));
            }

            if (!eta_coef_table.TryGetValue(0, out BesselYEtaTable eta_table)) {
                eta_table = new BesselYEtaTable(0);
                eta_coef_table.Add(0, eta_table);
            }

            BesselYEtaTable eta = eta_table;

            Complex<N> f0 = 1e-256, f1 = 0d, lambda = 0d;
            Complex<N> se = 0d;
            Complex<N> v = 1d / x;

            for (int k = m; k >= 1; k--) {
                if ((k & 1) == 0) {
                    lambda += f0;

                    se += f0 * eta[k / 2];
                }

                (f0, f1) = (2 * k * v * f0 - f1, f0);
            }

            lambda = Complex<N>.Ldexp(lambda, 1) + f0;

            Complex<N> y0 = Complex<N>.Ldexp((se + f0 * (Complex<N>.Log(Complex<N>.Ldexp(x, -1)) + MultiPrecision<N>.EulerGamma)) / (MultiPrecision<N>.PI * lambda), 1);

            return y0;
        }

        private static Complex<N> BesselY1Kernel(Complex<N> x, int m) {
            if (m < 2 || (m & 1) != 0) {
                throw new ArgumentOutOfRangeException(nameof(m));
            }

            if (!xi_coef_table.ContainsKey(0)) {
                if (!eta_coef_table.ContainsKey(0)) {
                    eta_coef_table.Add(0, new BesselYEtaTable(0));
                }

                xi_coef_table.Add(0, new BesselYXiTable(0, eta_coef_table[0]));
            }

            BesselYXiTable xi = xi_coef_table[0];

            Complex<N> f0 = 1e-256, f1 = 0d, lambda = 0d;
            Complex<N> sx = 0d;
            Complex<N> v = 1d / x;

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

            Complex<N> y1 = Complex<N>.Ldexp((sx - v * f0 + (Complex<N>.Log(Complex<N>.Ldexp(x, -1)) + MultiPrecision<N>.EulerGamma - 1d) * f1) / (lambda * MultiPrecision<N>.PI), 1);

            return y1;
        }

        internal static Complex<N> BesselIKernel(int n, Complex<N> x, int m, bool scale = false) {
            if (m < 2 || (m & 1) != 0 || n >= m) {
                throw new ArgumentOutOfRangeException(nameof(m));
            }

            n = int.Abs(n);

            if (n == 0) {
                return BesselI0Kernel(x, m, scale);
            }
            if (n == 1) {
                return BesselI1Kernel(x, m, scale);
            }

            Complex<N> f0 = 1e-256, f1 = 0d, lambda = 0d, fn = 0d;
            Complex<N> v = 1d / x;

            for (int k = m; k >= 1; k--) {
                lambda += f0;

                (f0, f1) = (2 * k * v * f0 + f1, f0);

                if (k - 1 == n) {
                    fn = f0;
                }
            }

            lambda = Complex<N>.Ldexp(lambda, 1) + f0;

            Complex<N> yn = fn / lambda;

            if (!scale) {
                yn *= Complex<N>.Exp(x);
            }

            return yn;
        }

        internal static Complex<N> BesselIKernel(MultiPrecision<N> nu, Complex<N> x, int m, bool scale = false) {
            int n = (int)MultiPrecision<N>.Floor(nu);
            MultiPrecision<N> alpha = nu - n;

            if (alpha == 0d) {
                return BesselIKernel(n, x, m, scale);
            }

            if (m < 2 || (m & 1) != 0 || n >= m) {
                throw new ArgumentOutOfRangeException(nameof(m));
            }

            if (!psi_coef_table.TryGetValue(alpha, out BesselIPsiTable psi_table)) {
                psi_table = new BesselIPsiTable(alpha);
                psi_coef_table.Add(alpha, psi_table);
            }

            BesselIPsiTable psi = psi_table;

            Complex<N> g0 = 1e-256, g1 = 0d, lambda = 0d;
            Complex<N> v = 1d / x;

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

                if (!scale) {
                    yn *= Complex<N>.Exp(x);
                }

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

                if (!scale) {
                    yn *= Complex<N>.Exp(x);
                }

                return yn;
            }
        }

        private static Complex<N> BesselI0Kernel(Complex<N> x, int m, bool scale = false) {
            if (m < 2 || (m & 1) != 0) {
                throw new ArgumentOutOfRangeException(nameof(m));
            }

            Complex<N> g0 = 1e-256, g1 = 0d, lambda = 0d;
            Complex<N> v = 1d / x;

            for (int k = m; k >= 1; k--) {
                lambda += g0;

                (g0, g1) = (2 * k * v * g0 + g1, g0);
            }

            lambda = Complex<N>.Ldexp(lambda, 1) + g0;

            Complex<N> y0 = g0 / lambda;

            if (!scale) {
                y0 *= Complex<N>.Exp(x);
            }

            return y0;
        }

        private static Complex<N> BesselI1Kernel(Complex<N> x, int m, bool scale = false) {
            if (m < 2 || (m & 1) != 0) {
                throw new ArgumentOutOfRangeException(nameof(m));
            }

            Complex<N> g0 = 1e-256, g1 = 0d, lambda = 0d;
            Complex<N> v = 1d / x;

            for (int k = m; k >= 1; k--) {
                lambda += g0;

                (g0, g1) = (2 * k * v * g0 + g1, g0);
            }

            lambda = Complex<N>.Ldexp(lambda, 1) + g0;

            Complex<N> y1 = g1 / lambda;

            if (!scale) {
                y1 *= Complex<N>.Exp(x);
            }

            return y1;
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

            public MultiPrecision<N> this[int n] => Value(n);

            private MultiPrecision<N> Value(int n) {
                ArgumentOutOfRangeException.ThrowIfNegative(n, nameof(n));

                if (n < table.Count) {
                    return table[n];
                }

                for (int m = table.Count; m <= n; m++) {
                    g = g * (alpha + m - 1d) / m;

                    MultiPrecision<N> phi = g * (alpha + 2 * m);

                    table.Add(phi);
                }

                return table[n];
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

            public MultiPrecision<N> this[int n] => Value(n);

            private MultiPrecision<N> Value(int n) {
                ArgumentOutOfRangeException.ThrowIfNegative(n, nameof(n));

                if (n < table.Count) {
                    return table[n];
                }

                for (int m = table.Count; m <= n; m++) {
                    g = g * (MultiPrecision<N>.Ldexp(alpha, 1) + m - 1d) / m;

                    MultiPrecision<N> phi = g * (alpha + m);

                    table.Add(phi);
                }

                return table[n];
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

            public MultiPrecision<N> this[int n] => Value(n);

            private MultiPrecision<N> Value(int n) {
                ArgumentOutOfRangeException.ThrowIfNegative(n, nameof(n));

                if (n < table.Count) {
                    return table[n];
                }

                for (int m = table.Count; m <= n; m++) {
                    if (alpha > 0d) {
                        g = -g * (alpha + m - 1) * (MultiPrecision<N>.Ldexp(alpha, 1) + m - 1d) / (m * (m - alpha));

                        MultiPrecision<N> eta = g * (alpha + 2 * m);

                        table.Add(eta);
                    }
                    else {
                        MultiPrecision<N> eta = (MultiPrecision<N>)2d / m;

                        table.Add(((m & 1) == 1) ? eta : -eta);
                    }
                }

                return table[n];
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

            public MultiPrecision<N> this[int n] => Value(n);

            private MultiPrecision<N> Value(int n) {
                ArgumentOutOfRangeException.ThrowIfNegative(n, nameof(n));

                if (n < table.Count) {
                    return table[n];
                }

                for (int m = table.Count; m <= n; m++) {
                    if (alpha > 0d) {
                        if ((m & 1) == 0) {
                            table.Add(eta[m / 2]);
                        }
                        else {
                            table.Add((eta[m / 2] - eta[m / 2 + 1]) / 2);
                        }
                    }
                    else {
                        if ((m & 1) == 1) {
                            MultiPrecision<N> xi = (MultiPrecision<N>)(2 * (m / 2) + 1) / (m / 2 * ((m / 2) + 1));
                            table.Add(((m & 2) > 0) ? xi : -xi);
                        }
                        else {
                            table.Add(MultiPrecision<N>.NaN);
                        }
                    }
                }

                return table[n];
            }
        };
    }
}
