﻿using MultiPrecision;
using System.Numerics;

namespace ComplexBessel {
    public static class YoshidaCoef<N> where N : struct, IConstant {
        struct Plus32 : IConstant {
            public readonly int Value => checked(default(N).Value + 32);
        }

        public static (MultiPrecision<N>[] cs, MultiPrecision<N>[] ds) Table(MultiPrecision<N> nu, int m) {
            Limit<Plus32>.HankelExpansion hankel = new(nu.Convert<Plus32>());

            MultiPrecision<N>[] cs = new MultiPrecision<N>[m + 1], ds = new MultiPrecision<N>[m + 1];

            for (int j = 0; j <= m; j++) {
                ds[j] = ShiftedLegendre.Table(m, m - j) / ((m - j + 1) * hankel.ACoef(m - j + 1)).Convert<N>();

                MultiPrecision<Plus32> sum = 0;

                for (int p = m - j; p <= m; p++) {
                    sum += ShiftedLegendre.Table(m, p) * hankel.ACoef(p - m + j) / ((p + 1) * hankel.ACoef(p + 1));
                }

                cs[j] = sum.Convert<N>();
            }

            return (cs, ds);
        }

        public static MultiPrecision<N>[][] Table(int m) {
            MultiPrecision<N>[][] dss = new MultiPrecision<N>[m + 1][];

            BigInteger[] fs = new BigInteger[m + 2];
            BigInteger[][] ps = new BigInteger[m + 1][], qs = new BigInteger[m + 1][];

            (fs[0], fs[1]) = (1, 1);
            (ps[0], ps[1]) = ([1], [-1, 1]);
            (qs[0], qs[1]) = ([1], [-(2 * m + 1) * (2 * m + 1), 1]);

            for (int k = 2; k <= m; k++) {
                (ps[k], qs[k]) = (new BigInteger[k + 1], new BigInteger[k + 1]);

                int sq2km1 = (2 * k - 1) * (2 * k - 1), sq2mkp3 = (2 * (m - k) + 3) * (2 * (m - k) + 3);

                ps[k][0] = -sq2km1 * ps[k - 1][0];
                ps[k][k] = 1;

                qs[k][0] = -sq2mkp3 * qs[k - 1][0];
                qs[k][k] = 1;

                for (int l = 1; l < k; l++) {
                    ps[k][l] = ps[k - 1][l - 1] - sq2km1 * ps[k - 1][l];
                    qs[k][l] = qs[k - 1][l - 1] - sq2mkp3 * qs[k - 1][l];
                }

                fs[k] = k * fs[k - 1];
            }

            fs[m + 1] = (m + 1) * fs[m];

            for (int i = 0; i <= m; i++) {
                dss[i] = new MultiPrecision<N>[i + 1];

                for (int j = 0; j <= i; j++) {
                    MultiPrecision<Plus32> b = 0;

                    for (int l = 0; l <= j; l++) {
                        for (int k = j - l; k <= i - l; k++) {
                            b += (MultiPrecision<Plus32>)(ShiftedLegendre.Table(m, m - k) * ps[i - k][l] * qs[k][j - l] * fs[m - k]) / (fs[i - k] * fs[m + 1]);
                        }
                    }

                    dss[i][j] = MultiPrecision<Plus32>.Ldexp(b, 2 * j - 3 * i).Convert<N>();
                }
            }

            return dss;
        }

        public static (MultiPrecision<N>[] cs, MultiPrecision<N>[] ds) Table(MultiPrecision<N> nu, MultiPrecision<N>[][] dss) {
            int m = dss.Length - 1;

            MultiPrecision<N> squa_nu = nu * nu;
            MultiPrecision<N>[] cs = new MultiPrecision<N>[m + 1], ds = new MultiPrecision<N>[m + 1];

            for (int i = 0; i <= m; i++) {
                MultiPrecision<N> d = dss[i][i], c = 0;

                for (int l = 0; l < i; l++) {
                    d *= MultiPrecision<N>.Square(m - l + MultiPrecision<N>.Point5) - squa_nu;
                }

                MultiPrecision<N> u = 1;
                for (int j = 0; j <= i; j++) {
                    c += dss[i][j] * u;
                    u *= squa_nu;
                }

                cs[i] = c;
                ds[i] = d;
            }

            return (cs, ds);
        }
    }
}
