using DoubleDouble;
using DoubleDoubleComplex;
using System.Collections.Concurrent;
using System.Collections.ObjectModel;
using System.Diagnostics;

namespace DDoubleComplexBessel {
    public class MillerBackward {
        public static readonly int BesselYEpsExponent = -12;

        private static readonly ConcurrentDictionary<ddouble, BesselJPhiTable> phi_coef_table = [];
        private static readonly ConcurrentDictionary<ddouble, BesselIPsiTable> psi_coef_table = [];
        private static readonly ConcurrentDictionary<ddouble, BesselYEtaTable> eta_coef_table = [];
        private static readonly ConcurrentDictionary<ddouble, BesselYXiTable> xi_coef_table = [];

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

        public static Complex BesselI(int n, Complex z) {
            Debug.Assert(ddouble.IsPositive(z.R));
            Debug.Assert(ddouble.IsPositive(z.I));

            int m = BesselIIterM((double)z.R, (double)z.I);

            Complex y = BesselIKernel(n, z, m);

            return y;
        }

        public static Complex BesselI(ddouble nu, Complex z) {
            Debug.Assert(ddouble.IsPositive(z.R));
            Debug.Assert(ddouble.IsPositive(z.I));

            int m = BesselIIterM((double)z.R, (double)z.I);

            if (BesselUtil.NearlyInteger(nu, out int n)) {
                Complex y = BesselIKernel(n, z, m);

                return y;
            }
            else {
                Complex y = BesselIKernel(nu, z, m);

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

        public static Complex BesselJKernel(ddouble nu, Complex z, int m) {
            int n = (int)ddouble.Floor(nu);
            ddouble alpha = nu - n;

            if (alpha == 0d) {
                return BesselJKernel(n, z, m);
            }

            Debug.Assert(m >= 2 && (m & 1) == 0 && n < m);

            if (!phi_coef_table.TryGetValue(alpha, out BesselJPhiTable phi)) {
                phi = new BesselJPhiTable(alpha);
                phi_coef_table[alpha] = phi;
            }

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

        private static Complex BesselJNKernel(int n, Complex z, int m) {
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

        private static Complex BesselYKernel(int n, Complex z, int m) {
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

        public static Complex BesselYKernel(ddouble nu, Complex z, int m) {
            int n = (int)ddouble.Round(nu);
            ddouble alpha = nu - n;

            if (alpha == 0d) {
                return BesselYKernel(n, z, m);
            }

            Debug.Assert(m >= 2 && (m & 1) == 0 && n < m);

            if (!eta_coef_table.TryGetValue(alpha, out BesselYEtaTable eta)) {
                eta = new BesselYEtaTable(alpha);
                eta_coef_table[alpha] = eta;
            }
            if (!xi_coef_table.TryGetValue(alpha, out BesselYXiTable xi)) {
                xi = new BesselYXiTable(alpha, eta);
                xi_coef_table[alpha] = xi;
            }
            if (!phi_coef_table.TryGetValue(alpha, out BesselJPhiTable phi)) {
                phi = new BesselJPhiTable(alpha);
                phi_coef_table[alpha] = phi;
            }

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

            ddouble rcot = 1d / ddouble.TanPi(alpha), rgamma = ddouble.Gamma(1d + alpha), rsqgamma = rgamma * rgamma;
            Complex r = Complex.Ldexp(ddouble.RcpPi * sqs, 1);
            Complex p = sqs * rsqgamma * ddouble.RcpPi;

            Complex xi0 = -Complex.Ldexp(v, 1) * p;

            (Complex eta0, Complex xi1) = ddouble.ILogB(alpha) >= BesselYEpsExponent
                ? (rcot - p / alpha, rcot + p * (alpha * (alpha + 1d) + 1d) / (alpha * (alpha - 1d)))
                : BesselYEta0Xi1Eps(alpha, z);

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

        private static (Complex eta0, Complex xi1) BesselYEta0Xi1Eps(ddouble alpha, Complex z) {
            const int N = 7;

            Complex lnz = Complex.Log(z);

            Complex eta0 = 0d, xi1 = 0d;
            for (int i = N, k = 1; i >= 0; i--) {
                Complex s = eta0_coef[^k], t = xi1_coef[^k];
                k++;

                for (int j = i; j >= 0; j--, k++) {
                    s = eta0_coef[k] + lnz * s;
                    t = xi1_coef[k] + lnz * t;
                }

                eta0 = s + alpha * eta0;
                xi1 = t + alpha * xi1;
            }

            return (eta0, xi1);
        }

        private static Complex BesselY0Kernel(Complex z, int m) {
            Debug.Assert(m >= 2 && (m & 1) == 0);

            if (!eta_coef_table.TryGetValue(0, out BesselYEtaTable eta)) {
                eta = new BesselYEtaTable(0);
                eta_coef_table[0] = eta;
            }

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

            Complex y0 = Complex.Ldexp((se + f0 * (Complex.Log(Complex.Ldexp(z, -1)) + ddouble.EulerGamma)) / (ddouble.Pi * lambda), 1);

            return y0;
        }

        private static Complex BesselY1Kernel(Complex z, int m) {
            Debug.Assert(m >= 2 && (m & 1) == 0);

            if (!xi_coef_table.TryGetValue(0, out BesselYXiTable xi)) {
                if (!eta_coef_table.TryGetValue(0, out BesselYEtaTable eta)) {
                    eta = new BesselYEtaTable(0);
                    eta_coef_table[0] = eta;
                }

                xi = new BesselYXiTable(0, eta);
                xi_coef_table[0] = xi;
            }

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

            Complex y1 = Complex.Ldexp((sx - v * f0 + (Complex.Log(Complex.Ldexp(z, -1)) + ddouble.EulerGamma - 1d) * f1) / (lambda * ddouble.Pi), 1);

            return y1;
        }

        private static Complex BesselYNKernel(int n, Complex z, int m) {
            if (!eta_coef_table.TryGetValue(0, out BesselYEtaTable eta)) {
                eta = new BesselYEtaTable(0);
                eta_coef_table[0] = eta;
            }

            if (!xi_coef_table.TryGetValue(0, out BesselYXiTable xi)) {
                xi = new BesselYXiTable(0, eta);
                xi_coef_table[0] = xi;
            }

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

            Complex yn = Complex.Ldexp(y1 / (lambda * ddouble.Pi), 1);

            return yn;
        }

        private static Complex BesselIKernel(int n, Complex z, int m) {
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

        public static Complex BesselIKernel(ddouble nu, Complex z, int m) {
            int n = (int)ddouble.Floor(nu);
            ddouble alpha = nu - n;

            if (alpha == 0d) {
                return BesselIKernel(n, z, m);
            }

            Debug.Assert(m >= 2 && (m & 1) == 0 && n < m);

            if (!psi_coef_table.TryGetValue(alpha, out BesselIPsiTable psi)) {
                psi = new BesselIPsiTable(alpha);
                psi_coef_table[alpha] = psi;
            }

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

                yn *= Complex.Exp(z);

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

                yn *= Complex.Exp(z);

                return yn;
            }
        }

        private static Complex BesselI0Kernel(Complex z, int m) {
            Debug.Assert(m >= 2 && (m & 1) == 0);

            Complex g0 = 1e-256, g1 = 0d, lambda = 0d;
            Complex v = 1d / z;

            for (int k = m; k >= 1; k--) {
                lambda += g0;

                (g0, g1) = (2 * k * v * g0 + g1, g0);
            }

            lambda = Complex.Ldexp(lambda, 1) + g0;

            Complex y0 = g0 / lambda;

            y0 *= Complex.Exp(z);

            return y0;
        }

        private static Complex BesselI1Kernel(Complex z, int m) {
            Debug.Assert(m >= 2 && (m & 1) == 0);

            Complex g0 = 1e-256, g1 = 0d, lambda = 0d;
            Complex v = 1d / z;

            for (int k = m; k >= 1; k--) {
                lambda += g0;

                (g0, g1) = (2 * k * v * g0 + g1, g0);
            }

            lambda = Complex.Ldexp(lambda, 1) + g0;

            Complex y1 = g1 / lambda;

            y1 *= Complex.Exp(z);

            return y1;
        }

        private static Complex BesselINKernel(int n, Complex z, int m) {
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

            yn *= Complex.Exp(z);

            return yn;
        }

        private class BesselJPhiTable {
            private readonly ddouble alpha;
            private readonly List<ddouble> table = [];

            private ddouble g;

            public BesselJPhiTable(ddouble alpha) {
                Debug.Assert(alpha > -1d && alpha < 1d, nameof(alpha));

                this.alpha = alpha;

                ddouble phi0 = ddouble.Gamma(1 + alpha);
                ddouble phi1 = phi0 * (alpha + 2d);

                this.g = phi0;

                this.table.Add(phi0);
                this.table.Add(phi1);
            }

            public ddouble this[int k] => Value(k);

            private ddouble Value(int k) {
                Debug.Assert(k >= 0);

                if (k < table.Count) {
                    return table[k];
                }

                lock (table) {
                    for (int i = table.Count; i <= k; i++) {
                        g = g * (alpha + i - 1d) / i;

                        ddouble phi = g * (alpha + 2 * i);

                        table.Add(phi);
                    }

                    return table[k];
                }
            }
        }

        private class BesselIPsiTable {
            private readonly ddouble alpha;
            private readonly List<ddouble> table = [];

            private ddouble g;

            public BesselIPsiTable(ddouble alpha) {
                Debug.Assert(alpha > -1d && alpha < 1d, nameof(alpha));

                this.alpha = alpha;

                ddouble psi0 = ddouble.Gamma(1d + alpha);
                ddouble psi1 = ddouble.Ldexp(psi0, 1) * (1d + alpha);

                this.g = ddouble.Ldexp(psi0, 1);

                this.table.Add(psi0);
                this.table.Add(psi1);
            }

            public ddouble this[int k] => Value(k);

            private ddouble Value(int k) {
                Debug.Assert(k >= 0);

                if (k < table.Count) {
                    return table[k];
                }

                lock (table) {
                    for (int i = table.Count; i <= k; i++) {
                        g = g * (ddouble.Ldexp(alpha, 1) + i - 1d) / i;

                        ddouble phi = g * (alpha + i);

                        table.Add(phi);
                    }

                    return table[k];
                }
            }
        }

        private class BesselYEtaTable {
            private readonly ddouble alpha;
            private readonly List<ddouble> table = [];

            private ddouble g;

            public BesselYEtaTable(ddouble alpha) {
                Debug.Assert(alpha > -1d && alpha < 1d, nameof(alpha));

                this.alpha = alpha;
                this.table.Add(ddouble.NaN);

                if (alpha != 0d) {
                    ddouble c = ddouble.Gamma(1d + alpha);
                    c *= c;
                    this.g = 1d / (1d - alpha) * c;

                    ddouble eta1 = (alpha + 2d) * g;

                    this.table.Add(eta1);
                }
            }

            public ddouble this[int k] => Value(k);

            private ddouble Value(int k) {
                Debug.Assert(k >= 0);

                if (k < table.Count) {
                    return table[k];
                }

                lock (table) {
                    for (int i = table.Count; i <= k; i++) {
                        if (alpha != 0d) {
                            g = -g * (alpha + i - 1) * (ddouble.Ldexp(alpha, 1) + i - 1d) / (i * (i - alpha));

                            ddouble eta = g * (alpha + 2 * i);

                            table.Add(eta);
                        }
                        else {
                            ddouble eta = (ddouble)2d / i;

                            table.Add((i & 1) == 1 ? eta : -eta);
                        }
                    }

                    return table[k];
                }
            }
        }

        private class BesselYXiTable {
            private readonly ddouble alpha;
            private readonly List<ddouble> table = [];
            private readonly BesselYEtaTable eta;

            public BesselYXiTable(ddouble alpha, BesselYEtaTable eta) {
                Debug.Assert(alpha >= -1d && alpha < 1d, nameof(alpha));

                this.alpha = alpha;
                this.table.Add(ddouble.NaN);
                this.table.Add(ddouble.NaN);

                this.eta = eta;
            }

            public ddouble this[int k] => Value(k);

            private ddouble Value(int k) {
                Debug.Assert(k >= 0);

                if (k < table.Count) {
                    return table[k];
                }

                lock (table) {
                    for (int i = table.Count; i <= k; i++) {
                        if (alpha != 0d) {
                            if ((i & 1) == 0) {
                                table.Add(eta[i / 2]);
                            }
                            else {
                                table.Add((eta[i / 2] - eta[i / 2 + 1]) / 2);
                            }
                        }
                        else {
                            if ((i & 1) == 1) {
                                ddouble xi = (ddouble)(2 * (i / 2) + 1) / (i / 2 * (i / 2 + 1));
                                table.Add((i & 2) > 0 ? xi : -xi);
                            }
                            else {
                                table.Add(ddouble.NaN);
                            }
                        }
                    }

                    return table[k];
                }
            }
        }

        private static readonly ReadOnlyCollection<ddouble> eta0_coef = new([
            (-1, -4, 0x9726B4CE5E80F444uL, 0x04F9CB1EFBE82ECCuL),
            (+1, -1, 0xA2F9836E4E441529uL, 0xFC2757D1F534DDC0uL),
            (-1, 0, 0xCA28399BC43B5702uL, 0xF0B44645EF33148BuL),
            (+1, -3, 0x9726B4CE5E80F444uL, 0x04F9CB1EFBE82ECCuL),
            (-1, -1, 0xA2F9836E4E441529uL, 0xFC2757D1F534DDC0uL),
            (+1, -3, 0x88365EC53E16C12CuL, 0x0368736407A27244uL),
            (+1, 0, 0x883B4FB4B14055BFuL, 0x85B55E7C87EB57FFuL),
            (-1, -3, 0x9726B4CE5E80F444uL, 0x04F9CB1EFBE82ECCuL),
            (+1, -2, 0xD94CAF3DBDB01C37uL, 0xFADF1FC29C467D01uL),
            (-1, 0, 0x9F9A4CA27EFFD464uL, 0xD17DF0A9838B4A3CuL),
            (-1, -2, 0x88365EC53E16C12CuL, 0x0368736407A27244uL),
            (-1, 0, 0x883B4FB4B14055BFuL, 0x85B55E7C87EB57FFuL),
            (+1, -4, 0xC988F11328ABF05AuL, 0xB14D0ED3FA8AE911uL),
            (-1, -3, 0xD94CAF3DBDB01C37uL, 0xFADF1FC29C467D01uL),
            (+1, -2, 0xD5CFA23DFF728E8CuL, 0x4FA05F9875BE7FDBuL),
            (+1, 0, 0x8ED06F76EF1F08D4uL, 0xFAAFF2FC72E0F22EuL),
            (+1, -2, 0x88365EC53E16C12CuL, 0x0368736407A27244uL),
            (+1, -1, 0xB5A46A4641AB1CFFuL, 0x5CF1D350B539CAA9uL),
            (-1, -5, 0xC988F11328ABF05AuL, 0xB14D0ED3FA8AE911uL),
            (+1, -4, 0xADD6F297CAF349C6uL, 0x624C19687D0530CDuL),
            (-1, 0, 0xA1EFBF6811197CC5uL, 0xEA0BD6E32998A1A1uL),
            (-1, -1, 0xD5CFA23DFF728E8CuL, 0x4FA05F9875BE7FDBuL),
            (-1, 0, 0x8ED06F76EF1F08D4uL, 0xFAAFF2FC72E0F22EuL),
            (-1, -3, 0xB59DD3B1A81E56E5uL, 0x59E099DAB4D8985AuL),
            (-1, -2, 0xB5A46A4641AB1CFFuL, 0x5CF1D350B539CAA9uL),
            (+1, -6, 0xA13A5A75BA2326AEuL, 0xF43DA5766208BA74uL),
            (-1, -6, 0xE7C943750E99B7B3uL, 0x2DBACC8B515C4112uL),
            (+1, -1, 0xA2EEA6FBB819E07EuL, 0x95CCA4D4A4C628CEuL),
            (+1, 0, 0x9E1267360C8533E2uL, 0xF1588EFC26ADAEDBuL),
            (+1, -1, 0xD5CFA23DFF728E8CuL, 0x4FA05F9875BE7FDBuL),
            (+1, -1, 0xBE6B3F493ED40BC6uL, 0xA39543FB43D6983EuL),
            (+1, -4, 0xB59DD3B1A81E56E5uL, 0x59E099DAB4D8985AuL),
            (+1, -3, 0x915055050155B0CCuL, 0x4A5B0F73C42E3BBAuL),
            (-1, -8, 0xD6F8789CF82EDE3EuL, 0x9AFCDC9DD80BA345uL),
            (+1, -7, 0x847301F9BF334466uL, 0x63462BBD5310252FuL),
            (-1, 0, 0xB0EDF45B2A0853AFuL, 0x2A42C979CA5D0596uL),
            (-1, 0, 0xA2EEA6FBB819E07EuL, 0x95CCA4D4A4C628CEuL),
            (-1, 0, 0x9E1267360C8533E2uL, 0xF1588EFC26ADAEDBuL),
            (-1, -1, 0x8E8A6C2954F709B2uL, 0xDFC03FBAF929AA92uL),
            (-1, -2, 0xBE6B3F493ED40BC6uL, 0xA39543FB43D6983EuL),
            (-1, -5, 0x914B0FC1534B78B7uL, 0x7B1A14AEF713AD15uL),
            (-1, -5, 0xC1C0715C01C79665uL, 0xB87969EFB03DA4F9uL),
            (+1, -10, 0xF5AE40B364C7D96CuL, 0x1ED7D78FD2567173uL),
            (-1, -9, 0x847301F9BF334466uL, 0x63462BBD5310252FuL),
        ]);

        private static readonly ReadOnlyCollection<ddouble> xi1_coef = new([
            (-1, -1, 0xB5DE5A081A1433B2uL, 0x7CC69135D4B1E39AuL),
            (+1, -1, 0xA2F9836E4E441529uL, 0xFC2757D1F534DDC0uL),
            (-1, 1, 0xABA41964255F42B5uL, 0x773880C3A34BE05AuL),
            (+1, 0, 0xB5DE5A081A1433B2uL, 0x7CC69135D4B1E39AuL),
            (-1, -1, 0xA2F9836E4E441529uL, 0xFC2757D1F534DDC0uL),
            (-1, 1, 0x86E3742ABAF45DA3uL, 0x21AA5401A70D1C66uL),
            (+1, 1, 0xD13DA106DF235947uL, 0xC0976A7F9B5A5829uL),
            (-1, 0, 0xB5DE5A081A1433B2uL, 0x7CC69135D4B1E39AuL),
            (+1, -2, 0xD94CAF3DBDB01C37uL, 0xFADF1FC29C467D01uL),
            (-1, 1, 0xF03C087CD2E2F132uL, 0xABA01CBF4A537060uL),
            (+1, 2, 0x86E3742ABAF45DA3uL, 0x21AA5401A70D1C66uL),
            (-1, 1, 0xD13DA106DF235947uL, 0xC0976A7F9B5A5829uL),
            (+1, -1, 0xF27DCD6022C59A43uL, 0x5108C19D1B97DA23uL),
            (-1, -3, 0xD94CAF3DBDB01C37uL, 0xFADF1FC29C467D01uL),
            (-1, 1, 0xC499BFB2F722CD86uL, 0x760E8ABF72CC4D3AuL),
            (+1, 2, 0xC422FE094F2AC935uL, 0x818D2129A54607CEuL),
            (-1, 2, 0x86E3742ABAF45DA3uL, 0x21AA5401A70D1C66uL),
            (+1, 1, 0x8B7E6B59EA1790DAuL, 0x8064F1AA6791901BuL),
            (-1, -2, 0xF27DCD6022C59A43uL, 0x5108C19D1B97DA23uL),
            (+1, -4, 0xADD6F297CAF349C6uL, 0x624C19687D0530CDuL),
            (-1, 2, 0x8F45E37E7DC47E26uL, 0x8FE6337E8ACC854DuL),
            (+1, 2, 0xC499BFB2F722CD86uL, 0x760E8ABF72CC4D3AuL),
            (-1, 2, 0xC422FE094F2AC935uL, 0x818D2129A54607CEuL),
            (+1, 1, 0xB3D9F038F945D22EuL, 0xD78DC5578966D089uL),
            (-1, 0, 0x8B7E6B59EA1790DAuL, 0x8064F1AA6791901BuL),
            (+1, -3, 0xC1FE3DE68237AE9CuL, 0x40D3CE174946481CuL),
            (-1, -6, 0xE7C943750E99B7B3uL, 0x2DBACC8B515C4112uL),
            (-1, 1, 0xD9277CDD4B4A0DEAuL, 0x740B8DE15C3A3708uL),
            (+1, 2, 0xF5188116761D8AE2uL, 0xE71C9F4A8A782581uL),
            (-1, 2, 0xC499BFB2F722CD86uL, 0x760E8ABF72CC4D3AuL),
            (+1, 2, 0x82C1FEB0DF71DB79uL, 0x0108C0C66E2EAFDEuL),
            (-1, 0, 0xB3D9F038F945D22EuL, 0xD78DC5578966D089uL),
            (+1, -2, 0xDF30ABC31025B490uL, 0xCD6E4F770C1C19C6uL),
            (-1, -4, 0x8154294456CFC9BDuL, 0x808D340F862EDABDuL),
            (+1, -7, 0x847301F9BF334466uL, 0x63462BBD5310252FuL),
            (-1, 2, 0x9833B38CBAB4864DuL, 0x9007F69410F14DABuL),
            (+1, 2, 0xD9277CDD4B4A0DEAuL, 0x740B8DE15C3A3708uL),
            (-1, 2, 0xF5188116761D8AE2uL, 0xE71C9F4A8A782581uL),
            (+1, 2, 0x83112A774F6C8904uL, 0x4EB45C7FA1DD88D1uL),
            (-1, 1, 0x82C1FEB0DF71DB79uL, 0x0108C0C66E2EAFDEuL),
            (+1, -1, 0x8FE18CFA6104A825uL, 0x793E37793AB8A6D4uL),
            (-1, -3, 0x94CB1D2CB56E7860uL, 0x88F434FA0812BBD9uL),
            (+1, -6, 0x93CDE604F57FC1FDuL, 0x2533A93650358C46uL),
            (-1, -9, 0x847301F9BF334466uL, 0x63462BBD5310252FuL),
        ]);
    }
}
