using DoubleDouble;
using MultiPrecision;
using MultiPrecisionComplex;
using System.Diagnostics;

namespace ComplexBesselSandbox {
    internal class Program {
        static void Main() {
            for (double alpha = -0.5; alpha <= 0.5; alpha += 1d / 64) {
                ddouble c0 = BesselYKernel(2 + alpha, 2);
                                
                Console.WriteLine($"{alpha},{c0}");
            }

            for (double alpha = -MillerBwdBesselYEps * 8; alpha <= MillerBwdBesselYEps * 8; alpha += MillerBwdBesselYEps / 8) {
                ddouble c0 = BesselYKernel(2 + alpha, 2);
                                
                Console.WriteLine($"{alpha},{c0}");
            }

            Console.WriteLine("END");
            Console.Read();
        }

        public static readonly double MillerBwdBesselYEps = double.ScaleB(1, -30);

        private static ddouble BesselYKernel(ddouble nu, ddouble x, int m = 256) {
            int n = (int)ddouble.Round(nu);
            ddouble alpha = nu - n;

            BesselYEtaTable eta = new BesselYEtaTable(alpha);
            BesselYXiTable xi = new BesselYXiTable(alpha, eta);
            BesselJPhiTable phi = new BesselJPhiTable(alpha);

            ddouble f0 = 1e-256, f1 = 0d, lambda = 0d;
            ddouble se = 0d, sxo = 0d, sxe = 0d;
            ddouble v = 1d / x;

            for (int k = m; k >= 1; k--) {
                if ((k & 1) == 0) {
                    lambda += f0 * phi[k / 2];

                    se += f0 * eta[k / 2];
                    sxe += f0 * xi[k];
                }
                else if (k >= 3) {
                    sxo += f0 * xi[k];
                }

                (f0, f1) = (ddouble.Ldexp(k + alpha, 1) * v * f0 - f1, f0);
            }

            ddouble s = ddouble.Pow(ddouble.Ldexp(v, 1), alpha), sqs = s * s;

            lambda += f0 * phi[0];
            lambda *= s;

            ddouble rcot = 1d / ddouble.TanPI(alpha), rgamma = ddouble.Gamma(1d + alpha), rsqgamma = rgamma * rgamma;
            ddouble r = ddouble.Ldexp(ddouble.RcpPI * sqs, 1);
            ddouble p = sqs * rsqgamma * ddouble.RcpPI;

            ddouble eta0 = ddouble.Abs(alpha) > MillerBwdBesselYEps
                ? rcot - p / alpha
                : BesselYEta0Eps(alpha, x);

            ddouble xi0 = -ddouble.Ldexp(v, 1) * p;
            ddouble xi1 = ddouble.Abs(alpha) > MillerBwdBesselYEps
                ? rcot + p * (alpha * (alpha + 1d) + 1d) / (alpha * (alpha - 1d))
                : BesselYXi1Eps(alpha, x);

            ddouble y0 = r * se + eta0 * f0;
            ddouble y1 = r * (3d * alpha * v * sxe + sxo) + xi0 * f0 + xi1 * f1;

            //
            //Console.WriteLine($"{nameof(n     )}: {n     }");
            //Console.WriteLine($"{nameof(alpha )}: {alpha }");
            //Console.WriteLine($"{nameof(f0    )}: {f0    }");
            //Console.WriteLine($"{nameof(f1    )}: {f1    }");
            //Console.WriteLine($"{nameof(eta0  )}: {eta0  }");
            //Console.WriteLine($"{nameof(xi0   )}: {xi0   }");
            //Console.WriteLine($"{nameof(xi1   )}: {xi1   }");
            //Console.WriteLine($"{nameof(y0    )}: {y0    }");
            //Console.WriteLine($"{nameof(y1    )}: {y1    }");
            //Console.WriteLine($"{nameof(lambda)}: {lambda}");
            //


            if (n == 0) {
                ddouble yn = y0 / lambda;

                return yn;
            }
            if (n == 1) {
                ddouble yn = y1 / lambda;

                return yn;
            }
            if (n >= 0) {
                for (int k = 1; k < n; k++) {
                    (y1, y0) = (ddouble.Ldexp(k + alpha, 1) * v * y1 - y0, y1);
                }

                ddouble yn = y1 / lambda;

                return yn;
            }
            else {
                for (int k = 0; k > n; k--) {
                    (y0, y1) = (ddouble.Ldexp(k + alpha, 1) * v * y0 - y1, y0);
                }

                ddouble yn = y0 / lambda;

                return yn;
            }
        }

        private static ddouble BesselYEta0Eps(ddouble alpha, ddouble x) {
            ddouble lnx = ddouble.Log(x), lnhalfx = ddouble.Log(ddouble.Ldexp(x, -1));
            ddouble pi = ddouble.PI, sqpi = pi * pi;
            ddouble ln2 = ddouble.Ln2, sqln2 = ln2 * ln2, cbln2 = sqln2 * ln2, qdln2 = sqln2 * sqln2;
            ddouble g = ddouble.EulerGamma;

            ddouble r0 = lnhalfx + g;
            ddouble r1 =
                (-sqln2 + lnx * (ln2 * 2d - lnx)) * 4d
                - sqpi
                - g * (lnhalfx * 2d + g) * 4d;
            ddouble r2 =
                (-cbln2 + lnx * (sqln2 * 3d + lnx * (ln2 * -3d + lnx))) * 4d
                + ddouble.Zeta3 * 2d
                + sqpi * (lnhalfx + g)
                + g * ((sqln2 + lnx * (ln2 * -2d + lnx)) * 3d + g * (lnhalfx * 3d + g)) * 4d;
            ddouble r3 =
                (-qdln2 + lnx * (cbln2 * 4d + lnx * (sqln2 * -6d + lnx * (ln2 * 4d - lnx)))) * 16d
                - ddouble.Zeta3 * (lnhalfx + g) * 32d
                - sqpi * ((sqln2 + lnx * (-ln2 * 2d + lnx) + g * (lnhalfx * 2d + g)) * 8d + sqpi)
                + g * ((cbln2 + lnx * (sqln2 * -3d + lnx * (ln2 * 3d - lnx))) * 4d
                + g * ((sqln2 + lnx * (ln2 * -2d + lnx)) * -6d
                + g * (lnhalfx * -4d
                - g))) * 16d;

            ddouble eta0 = (r0 * 48d + alpha * (r1 * 12d + alpha * (r2 * 8d + alpha * r3))) / (24d * ddouble.PI);

            return eta0;
        }

        private static ddouble BesselYXi1Eps(ddouble alpha, ddouble x) {
            ddouble lnx = ddouble.Log(x), lnhalfx = ddouble.Log(ddouble.Ldexp(x, -1)), lnxm1 = lnx - 1, lnhalfxm1 = lnhalfx - 1;
            ddouble pi = ddouble.PI, sqpi = pi * pi;
            ddouble ln2 = ddouble.Ln2, sqln2 = ln2 * ln2, cbln2 = sqln2 * ln2, qdln2 = sqln2 * sqln2;
            ddouble g = ddouble.EulerGamma;

            ddouble r0 = lnhalfxm1 + g;
            ddouble r1 =
                (-sqln2 + ln2 * lnxm1 * 2d + lnx * (2 - lnx)) * 4d
                - sqpi
                - g * (lnhalfxm1 * 2d + g) * 4d
                - 6d;
            ddouble r2 =
                -cbln2 * 4d + sqln2 * lnxm1 * 12d + lnx * (18d + lnx * (-12d + lnx * 4d))
                + ln2 * (lnx * (2d - lnx) * 12d - 18d)
                + ddouble.Zeta3 * 2d
                + sqpi * (lnhalfxm1 + g)
                + g * ((sqln2 - ln2 * lnxm1 * 2d + lnx * (-2d + lnx)) * 12d + 18d
                + g * (lnhalfxm1 * 12d
                + g * 4d))
                - 9d;
            ddouble r3 =
                -qdln2 * 16d
                + cbln2 * lnxm1 * 64d
                + sqln2 * (lnx * (2d - lnx) * 96d - 144d)
                + ln2 * (lnx * (9d + lnx * (-6d + lnx * 2d)) * 32d - 144d)
                + lnx * (9d + lnx * (-9d + lnx * (4d - lnx))) * 16d
                + ddouble.Zeta3 * (lnhalfxm1 + g) * -32d
                + sqpi * ((-sqln2 + ln2 * lnxm1 * 2d + lnx * (2d - lnx) - g * (lnhalfxm1 * 2d + g)) * 8d - 12d - sqpi)
                + g * ((cbln2 - sqln2 * lnxm1 * 3d) * 64d + ln2 * (lnx * (-2d + lnx) * 192d + 288d) + lnx * (-9d + lnx * (6d - lnx * 2d)) * 32d + 144d
                + g * ((-sqln2 + ln2 * lnxm1 * 2d + lnx * (2d - lnx)) * 96d - 144d
                + g * (lnhalfxm1 * -64d
                - g * 16d)))
                - 72d;

            ddouble xi1 = (r0 * 48d + alpha * (r1 * 12d + alpha * (r2 * 8d + alpha * r3))) / (24d * ddouble.PI);

            return xi1;
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

                for (int i = table.Count; i <= k; i++) {
                    g = g * (alpha + i - 1d) / i;

                    ddouble phi = g * (alpha + 2 * i);

                    table.Add(phi);
                }

                return table[k];
            }
        };

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

                for (int i = table.Count; i <= k; i++) {
                    g = g * (ddouble.Ldexp(alpha, 1) + i - 1d) / i;

                    ddouble phi = g * (alpha + i);

                    table.Add(phi);
                }

                return table[k];
            }
        };

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
        };

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
}
