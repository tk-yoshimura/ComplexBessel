using MultiPrecision;
using MultiPrecisionComplex;
using System.Diagnostics;

namespace ComplexBessel {
    public class MillerBackward<N> where N : struct, IConstant {
        struct DoubleN : IConstant {
            public readonly int Value => checked(default(N).Value * 2);
        }

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
            int n = (int)MultiPrecision<N>.Round(nu);
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

            MultiPrecision<DoubleN> rcot = 1d / MultiPrecision<DoubleN>.TanPI(alpha.Convert<DoubleN>());
            MultiPrecision<DoubleN> rgamma = MultiPrecision<DoubleN>.Gamma(1d + alpha.Convert<DoubleN>()), rsqgamma = rgamma * rgamma;
            Complex<N> r = Complex<DoubleN>.Ldexp(MultiPrecision<DoubleN>.RcpPI * sqs.Convert<DoubleN>(), 1).Convert<N>();
            Complex<DoubleN> p = sqs.Convert<DoubleN>() * rsqgamma * MultiPrecision<DoubleN>.RcpPI;

            Complex<N> xi0 = -Complex<N>.Ldexp(v, 1) * p.Convert<N>();

            Complex<N> eta0 = (rcot - p / alpha.Convert<DoubleN>()).Convert<N>();
            Complex<N> xi1 = (rcot + p * (alpha.Convert<DoubleN>() * (alpha.Convert<DoubleN>() + 1d) + 1d) / (alpha.Convert<DoubleN>() * (alpha.Convert<DoubleN>() - 1d))).Convert<N>();

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
                Debug.Assert(alpha > -1d && alpha < 1d, nameof(alpha));

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
        }

        private class BesselIPsiTable {
            private readonly MultiPrecision<N> alpha;
            private readonly List<MultiPrecision<N>> table = [];

            private MultiPrecision<N> g;

            public BesselIPsiTable(MultiPrecision<N> alpha) {
                Debug.Assert(alpha > -1d && alpha < 1d, nameof(alpha));

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
        }

        private class BesselYEtaTable {
            private readonly MultiPrecision<N> alpha;
            private readonly List<MultiPrecision<N>> table = [];

            private MultiPrecision<N> g;

            public BesselYEtaTable(MultiPrecision<N> alpha) {
                Debug.Assert(alpha > -1d && alpha < 1d, nameof(alpha));

                this.alpha = alpha;
                this.table.Add(MultiPrecision<N>.NaN);

                if (alpha != 0d) {
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
                    if (alpha != 0d) {
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
        }

        private class BesselYXiTable {
            private readonly MultiPrecision<N> alpha;
            private readonly List<MultiPrecision<N>> table = [];
            private readonly BesselYEtaTable eta;

            public BesselYXiTable(MultiPrecision<N> alpha, BesselYEtaTable eta) {
                Debug.Assert(alpha >= -1d && alpha < 1d, nameof(alpha));

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
