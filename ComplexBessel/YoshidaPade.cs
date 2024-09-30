﻿using MultiPrecision;
using MultiPrecisionComplex;
using System.Diagnostics;

namespace ComplexBessel {
    public class YoshidaPade<N> where N : struct, IConstant {
        const int m = 48;
        private static readonly MultiPrecision<N>[][] dss;

        public MultiPrecision<N> Nu { get; }
        private readonly MultiPrecision<N>[] cs0, ds0, cs1, ds1;

        static YoshidaPade() {
            dss = YoshidaCoef<N>.Table(m);
        }

        public YoshidaPade(MultiPrecision<N> nu) {
            Nu = MultiPrecision<N>.Abs(nu);

            (cs0, ds0) = YoshidaCoef<N>.Table(Nu - MultiPrecision<N>.Floor(Nu), dss);
            (cs1, ds1) = YoshidaCoef<N>.Table(Nu - MultiPrecision<N>.Floor(Nu) + 1, dss);
        }

        public Complex<N> BesselK(Complex<N> z) {
            Debug.Assert(z.R.Sign == Sign.Plus);
            Debug.Assert(z.I.Sign == Sign.Plus);

            if (Nu < 1d) {
                Complex<N> y = Value(z, cs0, ds0);

                return y;
            }
            else if (Nu < 2d) {
                Complex<N> y = Value(z, cs1, ds1);

                return y;
            }
            else {
                int n = (int)MultiPrecision<N>.Floor(Nu);
                MultiPrecision<N> alpha = Nu - n;

                Complex<N> y0 = Value(z, cs0, ds0);
                Complex<N> y1 = Value(z, cs1, ds1);

                Complex<N> v = 1d / z;

                for (int k = 1; k < n; k++) {
                    (y1, y0) = (Complex<N>.Ldexp(k + alpha, 1) * v * y1 + y0, y1);
                }

                return y1;
            }
        }

        static Complex<N> Value(Complex<N> z, MultiPrecision<N>[] cs, MultiPrecision<N>[] ds) {
            Complex<N> t = 1 / z, tn = 1;
            Complex<N> c = 0, d = 0;

            for (int j = 0; j < cs.Length; j++) {
                c += cs[j] * tn;
                d += ds[j] * tn;
                tn *= t;
            }

            Complex<N> y = Complex<N>.Sqrt(t * (MultiPrecision<N>.PI / 2)) * c / d;
            y *= Complex<N>.Exp(-z);

            return y;
        }
    }
}
