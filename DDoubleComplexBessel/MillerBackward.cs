using DoubleDouble;
using DoubleDoubleComplex;
using System.Diagnostics;

namespace DDoubleComplexBessel {
    public class MillerBackward {
        private static readonly Dictionary<ddouble, BesselJPhiTable> phi_coef_table = [];
        private static readonly Dictionary<ddouble, BesselIPsiTable> psi_coef_table = [];
        private static readonly Dictionary<ddouble, BesselYEtaTable> eta_coef_table = [];
        private static readonly Dictionary<ddouble, BesselYXiTable> xi_coef_table = [];

        public static Complex BesselJ(int n, Complex z) {
            Debug.Assert(ddouble.IsPositive(z.R));
            Debug.Assert(ddouble.IsPositive(z.I));

            int m = BesselJYIterM((double)z.R);

            Complex y = BesselJKernel(n, z, m);

            return y;
        }

        public static Complex BesselJ(ddouble nu, Complex z) {
            Debug.Assert(ddouble.IsPositive(z.R));
            Debug.Assert(ddouble.IsPositive(z.I));

            int m = BesselJYIterM((double)z.R);

            if (BesselUtil.NearlyInteger(nu, out int n)) {
                Complex y = BesselJKernel(n, z, m);

                return y;
            }
            else {
                Complex y = BesselJKernel(nu, z, m);

                return y;
            }
        }

        public static Complex BesselY(int n, Complex z) {
            Debug.Assert(ddouble.IsPositive(z.R));
            Debug.Assert(ddouble.IsPositive(z.I));

            int m = BesselJYIterM((double)z.R);

            Complex y = BesselYKernel(n, z, m);

            return y;
        }

        public static Complex BesselY(ddouble nu, Complex z) {
            Debug.Assert(ddouble.IsPositive(z.R));
            Debug.Assert(ddouble.IsPositive(z.I));

            int m = BesselJYIterM((double)z.R);

            if (BesselUtil.NearlyInteger(nu, out int n)) {
                Complex y = BesselYKernel(n, z, m);

                return y;
            }
            else {
                Complex y = BesselYKernel(nu, z, m);

                return y;
            }
        }

        private static int BesselJYIterM(double r) {
            int m = (int)double.Ceiling(3.8029e1 + r * 1.6342e0);

            return (m + 1) / 2 * 2;
        }

        public static Complex BesselI(int n, Complex z, bool scale = false) {
            Debug.Assert(ddouble.IsPositive(z.R));
            Debug.Assert(ddouble.IsPositive(z.I));

            int m = BesselIIterM((double)z.R, (double)z.I);

            Complex y = BesselIKernel(n, z, m, scale);

            return y;
        }

        public static Complex BesselI(ddouble nu, Complex z, bool scale = false) {
            Debug.Assert(ddouble.IsPositive(z.R));
            Debug.Assert(ddouble.IsPositive(z.I));

            int m = BesselIIterM((double)z.R, (double)z.I);

            if (BesselUtil.NearlyInteger(nu, out int n)) {
                Complex y = BesselIKernel(n, z, m, scale);

                return y;
            }
            else {
                Complex y = BesselIKernel(nu, z, m, scale);

                return y;
            }
        }

        private static int BesselIIterM(double r, double i) {
            int m = (int)double.Ceiling(3.3612e1 + r * 1.3557e0 + i * 1.8485e0 - r * i * 4.3649e-2);

            return (m + 1) / 2 * 2;
        }

        private static Complex BesselJKernel(int n, Complex z, int m) {
            Debug.Assert(m >= 2 && (m & 1) == 0 && n < m);

            if (n < 0) {
                return ((n & 1) == 0) ? BesselJKernel(-n, z, m) : -BesselJKernel(-n, z, m);
            }
            if (n == 0) {
                return BesselJ0Kernel(z, m);
            }
            if (n == 1) {
                return BesselJ1Kernel(z, m);
            }

            Complex f0 = 1e-256d, f1 = 0d, fn = 0d, lambda = 0d;
            Complex v = 1d / z;

            for (int k = m; k >= 1; k--) {
                if ((k & 1) == 0) {
                    lambda += f0;
                }

                (f0, f1) = (2 * k * v * f0 - f1, f0);

                if (k - 1 == n) {
                    fn = f0;
                }
            }

            lambda = Complex.Ldexp(lambda, 1) + f0;

            Complex yn = fn / lambda;

            return yn;
        }

        public static Complex BesselJKernel(ddouble nu, Complex z, int m) {
            int n = (int)ddouble.Floor(nu);
            ddouble alpha = nu - n;

            if (alpha == 0d) {
                return BesselJKernel(n, z, m);
            }

            Debug.Assert(m >= 2 && (m & 1) == 0 && n < m);

            if (!phi_coef_table.TryGetValue(alpha, out BesselJPhiTable phi_table)) {
                phi_table = new BesselJPhiTable(alpha);
                phi_coef_table.Add(alpha, phi_table);
            }

            BesselJPhiTable phi = phi_table;

            Complex f0 = 1e-256d, f1 = 0d, lambda = 0d;
            Complex v = 1d / z;

            if (n >= 0) {
                Complex fn = 0d;

                for (int k = m; k >= 1; k--) {
                    if ((k & 1) == 0) {
                        lambda += f0 * phi[k / 2];
                    }

                    (f0, f1) = (Complex.Ldexp(k + alpha, 1) * v * f0 - f1, f0);

                    if (k - 1 == n) {
                        fn = f0;
                    }
                }

                lambda += f0 * phi[0];
                lambda *= Complex.Pow(Complex.Ldexp(v, 1), alpha);

                Complex yn = fn / lambda;

                return yn;
            }
            else {
                for (int k = m; k >= 1; k--) {
                    if ((k & 1) == 0) {
                        lambda += f0 * phi[k / 2];
                    }

                    (f0, f1) = (Complex.Ldexp(k + alpha, 1) * v * f0 - f1, f0);
                }

                lambda += f0 * phi[0];
                lambda *= Complex.Pow(Complex.Ldexp(v, 1), alpha);

                for (int k = 0; k > n; k--) {
                    (f0, f1) = (Complex.Ldexp(k + alpha, 1) * v * f0 - f1, f0);
                }

                Complex yn = f0 / lambda;

                return yn;
            }
        }

        private static Complex BesselJ0Kernel(Complex z, int m) {
            Debug.Assert(m >= 2 && (m & 1) == 0);

            Complex f0 = 1e-256d, f1 = 0d, lambda = 0d;
            Complex v = 1d / z;

            for (int k = m; k >= 1; k--) {
                if ((k & 1) == 0) {
                    lambda += f0;
                }

                (f0, f1) = (2 * k * v * f0 - f1, f0);
            }

            lambda = Complex.Ldexp(lambda, 1) + f0;

            Complex y0 = f0 / lambda;

            return y0;
        }

        private static Complex BesselJ1Kernel(Complex z, int m) {
            Debug.Assert(m >= 2 && (m & 1) == 0);

            Complex f0 = 1e-256d, f1 = 0d, lambda = 0d;
            Complex v = 1d / z;

            for (int k = m; k >= 1; k--) {
                if ((k & 1) == 0) {
                    lambda += f0;
                }

                (f0, f1) = (2 * k * v * f0 - f1, f0);
            }

            lambda = Complex.Ldexp(lambda, 1) + f0;

            Complex y1 = f1 / lambda;

            return y1;
        }

        private static Complex BesselYKernel(int n, Complex z, int m) {
            Debug.Assert(m >= 2 && (m & 1) == 0 && n < m);

            if (n < 0) {
                return ((n & 1) == 0) ? BesselYKernel(-n, z, m) : -BesselYKernel(-n, z, m);
            }
            if (n == 0) {
                return BesselY0Kernel(z, m);
            }
            if (n == 1) {
                return BesselY1Kernel(z, m);
            }

            if (!eta_coef_table.ContainsKey(0)) {
                eta_coef_table.Add(0, new BesselYEtaTable(0));
            }

            BesselYEtaTable eta = eta_coef_table[0];

            if (!xi_coef_table.ContainsKey(0)) {
                xi_coef_table.Add(0, new BesselYXiTable(0, eta));
            }

            BesselYXiTable xi = xi_coef_table[0];

            Complex f0 = 1e-256, f1 = 0d, lambda = 0d;
            Complex se = 0d, sx = 0d;
            Complex v = 1d / z;

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

            lambda = Complex.Ldexp(lambda, 1) + f0;

            Complex c = Complex.Log(Complex.Ldexp(z, -1)) + ddouble.EulerGamma;

            Complex y0 = se + f0 * c;
            Complex y1 = sx - v * f0 + (c - 1d) * f1;

            for (int k = 1; k < n; k++) {
                (y1, y0) = (2 * k * v * y1 - y0, y1);
            }

            Complex yn = Complex.Ldexp(y1 / (lambda * ddouble.PI), 1);

            return yn;
        }

        public static Complex BesselYKernel(ddouble nu, Complex z, int m) {
            int n = (int)ddouble.Floor(nu);
            ddouble alpha = nu - n;

            if (alpha == 0d) {
                return BesselYKernel(n, z, m);
            }

            Debug.Assert(m >= 2 && (m & 1) == 0 && n < m);

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

            Complex f0 = 1e-256, f1 = 0d, lambda = 0d;
            Complex se = 0d, sxo = 0d, sxe = 0d;
            Complex v = 1d / z;

            for (int k = m; k >= 1; k--) {
                if ((k & 1) == 0) {
                    lambda += f0 * phi[k / 2];

                    se += f0 * eta[k / 2];
                    sxe += f0 * xi[k];
                }
                else if (k >= 3) {
                    sxo += f0 * xi[k];
                }

                (f0, f1) = (Complex.Ldexp(k + alpha, 1) * v * f0 - f1, f0);
            }

            Complex s = Complex.Pow(Complex.Ldexp(v, 1), alpha), sqs = s * s;

            lambda += f0 * phi[0];
            lambda *= s;

            ddouble rcot = 1d / ddouble.TanPI(alpha), rgamma = ddouble.Gamma(1d + alpha), rsqgamma = rgamma * rgamma;
            Complex r = Complex.Ldexp(ddouble.RcpPI * sqs, 1);
            Complex p = sqs * rsqgamma * ddouble.RcpPI;

            Complex eta0 = (ddouble.Abs(alpha) > BesselUtil.MillerBwdBesselYEps)
                ? (rcot - p / alpha)
                : BesselYEta0Eps(alpha, z);

            Complex xi0 = -Complex.Ldexp(v, 1) * p;
            Complex xi1 = (ddouble.Abs(alpha) > BesselUtil.MillerBwdBesselYEps)
                ? rcot + p * (alpha * (alpha + 1d) + 1d) / (alpha * (alpha - 1d))
                : BesselYXi1Eps(alpha, z);

            Complex y0 = r * se + eta0 * f0;
            Complex y1 = r * (3d * alpha * v * sxe + sxo) + xi0 * f0 + xi1 * f1;

            if (n == 0) {
                Complex yn = y0 / lambda;

                return yn;
            }
            if (n == 1) {
                Complex yn = y1 / lambda;

                return yn;
            }
            if (n >= 0) {
                for (int k = 1; k < n; k++) {
                    (y1, y0) = (Complex.Ldexp(k + alpha, 1) * v * y1 - y0, y1);
                }

                Complex yn = y1 / lambda;

                return yn;
            }
            else {
                for (int k = 0; k > n; k--) {
                    (y0, y1) = (Complex.Ldexp(k + alpha, 1) * v * y0 - y1, y0);
                }

                Complex yn = y0 / lambda;

                return yn;
            }
        }

        private static Complex BesselYEta0Eps(ddouble alpha, Complex z) {
            Complex lnz = Complex.Log(z), lnhalfz = Complex.Log(Complex.Ldexp(z, -1));
            ddouble pi = ddouble.PI, sqpi = pi * pi;
            ddouble ln2 = ddouble.Ln2, sqln2 = ln2 * ln2, cbln2 = sqln2 * ln2, qdln2 = sqln2 * sqln2;
            ddouble g = ddouble.EulerGamma;

            Complex r0 = lnhalfz + g;
            Complex r1 =
                (-sqln2 + lnz * (ln2 * 2d - lnz)) * 4d
                - sqpi
                - g * (lnhalfz * 2d + g) * 4d;
            Complex r2 =
                (-cbln2 + lnz * (sqln2 * 3d + lnz * (ln2 * -3d + lnz))) * 4d
                + ddouble.Zeta3 * 2d
                + sqpi * (lnhalfz + g)
                + g * ((sqln2 + lnz * (ln2 * -2d + lnz)) * 3d + g * (lnhalfz * 3d + g)) * 4d;
            Complex r3 =
                (-qdln2 + lnz * (cbln2 * 4d + lnz * (sqln2 * -6d + lnz * (ln2 * 4d - lnz)))) * 16d
                - ddouble.Zeta3 * (lnhalfz + g) * 32d
                - sqpi * (((sqln2 + lnz * (-ln2 * 2d + lnz) + g * (lnhalfz * 2d + g)) * 8d) + sqpi)
                + g * ((cbln2 + lnz * (sqln2 * -3d + lnz * (ln2 * 3d - lnz))) * 4d
                + g * ((sqln2 + lnz * (ln2 * -2d + lnz)) * -6d
                + g * (lnhalfz * -4d
                - g))) * 16d;

            Complex eta0 = (r0 * 48d + alpha * (r1 * 12d + alpha * (r2 * 8d + alpha * r3))) / (24d * ddouble.PI);

            return eta0;
        }

        static Complex BesselYXi1Eps(ddouble alpha, Complex z) {
            Complex lnz = Complex.Log(z), lnhalfz = Complex.Log(Complex.Ldexp(z, -1)), lnxm1 = lnz - 1, lnhalfxm1 = lnhalfz - 1;
            ddouble pi = ddouble.PI, sqpi = pi * pi;
            ddouble ln2 = ddouble.Ln2, sqln2 = ln2 * ln2, cbln2 = sqln2 * ln2, qdln2 = sqln2 * sqln2;
            ddouble g = ddouble.EulerGamma;

            Complex r0 = lnhalfxm1 + g;
            Complex r1 =
                (-sqln2 + ln2 * lnxm1 * 2d + lnz * (2 - lnz)) * 4d
                - sqpi
                - g * (lnhalfxm1 * 2d + g) * 4d
                - 6d;
            Complex r2 =
                -cbln2 * 4d + sqln2 * lnxm1 * 12d + lnz * (18d + lnz * (-12d + lnz * 4d))
                + ln2 * (lnz * (2d - lnz) * 12d - 18d)
                + ddouble.Zeta3 * 2d
                + sqpi * (lnhalfxm1 + g)
                + g * ((sqln2 - ln2 * lnxm1 * 2d + lnz * (-2d + lnz)) * 12d + 18d
                + g * (lnhalfxm1 * 12d
                + g * 4d))
                - 9d;
            Complex r3 =
                -qdln2 * 16d
                + cbln2 * lnxm1 * 64d
                + sqln2 * (lnz * (2d - lnz) * 96d - 144d)
                + ln2 * (lnz * (9d + lnz * (-6d + lnz * 2d)) * 32d - 144d)
                + lnz * (9d + lnz * (-9d + lnz * (4d - lnz))) * 16d
                + ddouble.Zeta3 * (lnhalfxm1 + g) * -32d
                + sqpi * ((-sqln2 + ln2 * lnxm1 * 2d + lnz * (2d - lnz) - g * (lnhalfxm1 * 2d + g)) * 8d - 12d - sqpi)
                + g * ((cbln2 - sqln2 * lnxm1 * 3d) * 64d + ln2 * (lnz * (-2d + lnz) * 192d + 288d) + lnz * (-9d + lnz * (6d - lnz * 2d)) * 32d + 144d
                + g * ((-sqln2 + ln2 * lnxm1 * 2d + lnz * (2d - lnz)) * 96d - 144d
                + g * (lnhalfxm1 * -64d
                - g * 16d)))
                - 72d;

            Complex xi1 = (r0 * 48d + alpha * (r1 * 12d + alpha * (r2 * 8d + alpha * r3))) / (24d * ddouble.PI);

            return xi1;
        }

        private static Complex BesselY0Kernel(Complex z, int m) {
            Debug.Assert(m >= 2 && (m & 1) == 0);

            if (!eta_coef_table.TryGetValue(0, out BesselYEtaTable eta_table)) {
                eta_table = new BesselYEtaTable(0);
                eta_coef_table.Add(0, eta_table);
            }

            BesselYEtaTable eta = eta_table;

            Complex f0 = 1e-256, f1 = 0d, lambda = 0d;
            Complex se = 0d;
            Complex v = 1d / z;

            for (int k = m; k >= 1; k--) {
                if ((k & 1) == 0) {
                    lambda += f0;

                    se += f0 * eta[k / 2];
                }

                (f0, f1) = (2 * k * v * f0 - f1, f0);
            }

            lambda = Complex.Ldexp(lambda, 1) + f0;

            Complex y0 = Complex.Ldexp((se + f0 * (Complex.Log(Complex.Ldexp(z, -1)) + ddouble.EulerGamma)) / (ddouble.PI * lambda), 1);

            return y0;
        }

        private static Complex BesselY1Kernel(Complex z, int m) {
            Debug.Assert(m >= 2 && (m & 1) == 0);

            if (!xi_coef_table.ContainsKey(0)) {
                if (!eta_coef_table.ContainsKey(0)) {
                    eta_coef_table.Add(0, new BesselYEtaTable(0));
                }

                xi_coef_table.Add(0, new BesselYXiTable(0, eta_coef_table[0]));
            }

            BesselYXiTable xi = xi_coef_table[0];

            Complex f0 = 1e-256, f1 = 0d, lambda = 0d;
            Complex sx = 0d;
            Complex v = 1d / z;

            for (int k = m; k >= 1; k--) {
                if ((k & 1) == 0) {
                    lambda += f0;
                }
                else if (k >= 3) {
                    sx += f0 * xi[k];
                }

                (f0, f1) = (2 * k * v * f0 - f1, f0);
            }

            lambda = Complex.Ldexp(lambda, 1) + f0;

            Complex y1 = Complex.Ldexp((sx - v * f0 + (Complex.Log(Complex.Ldexp(z, -1)) + ddouble.EulerGamma - 1d) * f1) / (lambda * ddouble.PI), 1);

            return y1;
        }

        private static Complex BesselIKernel(int n, Complex z, int m, bool scale = false) {
            Debug.Assert(m >= 2 && (m & 1) == 0 && n < m);

            n = int.Abs(n);

            if (n == 0) {
                return BesselI0Kernel(z, m, scale);
            }
            if (n == 1) {
                return BesselI1Kernel(z, m, scale);
            }

            Complex f0 = 1e-256, f1 = 0d, lambda = 0d, fn = 0d;
            Complex v = 1d / z;

            for (int k = m; k >= 1; k--) {
                lambda += f0;

                (f0, f1) = (2 * k * v * f0 + f1, f0);

                if (k - 1 == n) {
                    fn = f0;
                }
            }

            lambda = Complex.Ldexp(lambda, 1) + f0;

            Complex yn = fn / lambda;

            if (!scale) {
                yn *= Complex.Exp(z);
            }

            return yn;
        }

        public static Complex BesselIKernel(ddouble nu, Complex z, int m, bool scale = false) {
            int n = (int)ddouble.Floor(nu);
            ddouble alpha = nu - n;

            if (alpha == 0d) {
                return BesselIKernel(n, z, m, scale);
            }

            Debug.Assert(m >= 2 && (m & 1) == 0 && n < m);

            if (!psi_coef_table.TryGetValue(alpha, out BesselIPsiTable psi_table)) {
                psi_table = new BesselIPsiTable(alpha);
                psi_coef_table.Add(alpha, psi_table);
            }

            BesselIPsiTable psi = psi_table;

            Complex g0 = 1e-256, g1 = 0d, lambda = 0d;
            Complex v = 1d / z;

            if (n >= 0) {
                Complex gn = 0d;

                for (int k = m; k >= 1; k--) {
                    lambda += g0 * psi[k];

                    (g0, g1) = (Complex.Ldexp(k + alpha, 1) * v * g0 + g1, g0);

                    if (k - 1 == n) {
                        gn = g0;
                    }
                }

                lambda += g0 * psi[0];
                lambda *= Complex.Pow(Complex.Ldexp(v, 1), alpha);

                Complex yn = gn / lambda;

                if (!scale) {
                    yn *= Complex.Exp(z);
                }

                return yn;
            }
            else {
                for (int k = m; k >= 1; k--) {
                    lambda += g0 * psi[k];

                    (g0, g1) = (Complex.Ldexp(k + alpha, 1) * v * g0 + g1, g0);
                }

                lambda += g0 * psi[0];
                lambda *= Complex.Pow(Complex.Ldexp(v, 1), alpha);

                for (int k = 0; k > n; k--) {
                    (g0, g1) = (Complex.Ldexp(k + alpha, 1) * v * g0 + g1, g0);
                }

                Complex yn = g0 / lambda;

                if (!scale) {
                    yn *= Complex.Exp(z);
                }

                return yn;
            }
        }

        private static Complex BesselI0Kernel(Complex z, int m, bool scale = false) {
            Debug.Assert(m >= 2 && (m & 1) == 0);

            Complex g0 = 1e-256, g1 = 0d, lambda = 0d;
            Complex v = 1d / z;

            for (int k = m; k >= 1; k--) {
                lambda += g0;

                (g0, g1) = (2 * k * v * g0 + g1, g0);
            }

            lambda = Complex.Ldexp(lambda, 1) + g0;

            Complex y0 = g0 / lambda;

            if (!scale) {
                y0 *= Complex.Exp(z);
            }

            return y0;
        }

        private static Complex BesselI1Kernel(Complex z, int m, bool scale = false) {
            Debug.Assert(m >= 2 && (m & 1) == 0);

            Complex g0 = 1e-256, g1 = 0d, lambda = 0d;
            Complex v = 1d / z;

            for (int k = m; k >= 1; k--) {
                lambda += g0;

                (g0, g1) = (2 * k * v * g0 + g1, g0);
            }

            lambda = Complex.Ldexp(lambda, 1) + g0;

            Complex y1 = g1 / lambda;

            if (!scale) {
                y1 *= Complex.Exp(z);
            }

            return y1;
        }

        private class BesselJPhiTable {
            private readonly ddouble alpha;
            private readonly List<ddouble> table = [];

            private ddouble g;

            public BesselJPhiTable(ddouble alpha) {
                Debug.Assert(alpha > 0d && alpha < 1d, nameof(alpha));

                this.alpha = alpha;

                ddouble phi0 = ddouble.Gamma(1 + alpha);
                ddouble phi1 = phi0 * (alpha + 2d);

                this.g = phi0;

                this.table.Add(phi0);
                this.table.Add(phi1);
            }

            public ddouble this[int n] => Value(n);

            private ddouble Value(int n) {
                Debug.Assert(n >= 0);

                if (n < table.Count) {
                    return table[n];
                }

                for (int m = table.Count; m <= n; m++) {
                    g = g * (alpha + m - 1d) / m;

                    ddouble phi = g * (alpha + 2 * m);

                    table.Add(phi);
                }

                return table[n];
            }
        };

        private class BesselIPsiTable {
            private readonly ddouble alpha;
            private readonly List<ddouble> table = [];

            private ddouble g;

            public BesselIPsiTable(ddouble alpha) {
                Debug.Assert(alpha > 0d && alpha < 1d, nameof(alpha));

                this.alpha = alpha;

                ddouble psi0 = ddouble.Gamma(1d + alpha);
                ddouble psi1 = ddouble.Ldexp(psi0, 1) * (1d + alpha);

                this.g = ddouble.Ldexp(psi0, 1);

                this.table.Add(psi0);
                this.table.Add(psi1);
            }

            public ddouble this[int n] => Value(n);

            private ddouble Value(int n) {
                Debug.Assert(n >= 0);

                if (n < table.Count) {
                    return table[n];
                }

                for (int m = table.Count; m <= n; m++) {
                    g = g * (ddouble.Ldexp(alpha, 1) + m - 1d) / m;

                    ddouble phi = g * (alpha + m);

                    table.Add(phi);
                }

                return table[n];
            }
        };

        private class BesselYEtaTable {
            private readonly ddouble alpha;
            private readonly List<ddouble> table = [];

            private ddouble g;

            public BesselYEtaTable(ddouble alpha) {
                Debug.Assert(alpha >= 0d && alpha < 1d, nameof(alpha));

                this.alpha = alpha;
                this.table.Add(ddouble.NaN);

                if (alpha > 0d) {
                    ddouble c = ddouble.Gamma(1d + alpha);
                    c *= c;
                    this.g = 1d / (1d - alpha) * c;

                    ddouble eta1 = (alpha + 2d) * g;

                    this.table.Add(eta1);
                }
            }

            public ddouble this[int n] => Value(n);

            private ddouble Value(int n) {
                Debug.Assert(n >= 0);

                if (n < table.Count) {
                    return table[n];
                }

                for (int m = table.Count; m <= n; m++) {
                    if (alpha > 0d) {
                        g = -g * (alpha + m - 1) * (ddouble.Ldexp(alpha, 1) + m - 1d) / (m * (m - alpha));

                        ddouble eta = g * (alpha + 2 * m);

                        table.Add(eta);
                    }
                    else {
                        ddouble eta = (ddouble)2d / m;

                        table.Add(((m & 1) == 1) ? eta : -eta);
                    }
                }

                return table[n];
            }
        };

        private class BesselYXiTable {
            private readonly ddouble alpha;
            private readonly List<ddouble> table = [];
            private readonly BesselYEtaTable eta;

            public BesselYXiTable(ddouble alpha, BesselYEtaTable eta) {
                Debug.Assert(alpha >= 0d && alpha < 1d, nameof(alpha));

                this.alpha = alpha;
                this.table.Add(ddouble.NaN);
                this.table.Add(ddouble.NaN);

                this.eta = eta;
            }

            public ddouble this[int n] => Value(n);

            private ddouble Value(int n) {
                Debug.Assert(n >= 0);

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
                            ddouble xi = (ddouble)(2 * (m / 2) + 1) / (m / 2 * ((m / 2) + 1));
                            table.Add(((m & 2) > 0) ? xi : -xi);
                        }
                        else {
                            table.Add(ddouble.NaN);
                        }
                    }
                }

                return table[n];
            }
        }
    }
}
