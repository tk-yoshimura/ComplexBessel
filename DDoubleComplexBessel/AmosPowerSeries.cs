﻿using DoubleDouble;
using DoubleDoubleComplex;
using System.Collections.ObjectModel;
using System.Diagnostics;

namespace DDoubleComplexBessel {
    public static class AmosPowerSeries {

        public static Complex BesselK(ddouble nu, Complex z) {
            Debug.Assert(nu >= 0d);

            int n = (int)ddouble.Round(nu);
            ddouble alpha = nu - n;

            Complex k0 = BesselKNearZeroNu(alpha, z, terms: 256);

            if (n == 0) {
                return k0;
            }

            Complex i0 = PowerSeries.BesselI(alpha, z), i1 = PowerSeries.BesselI(alpha + 1d, z);

            Complex k1 = (1d - i1 * k0 * z) / (i0 * z);

            if (n == 1) {
                return k1;
            }

            Complex v = 1d / z;

            for (int k = 1; k < n; k++) {
                (k1, k0) = (ddouble.Ldexp(k + alpha, 1) * v * k1 + k0, k1);
            }

            return k1;
        }

        private static Complex BesselKNearZeroNu(ddouble nu, Complex z, int terms) {
            nu = ddouble.Abs(nu);

            ddouble nu2 = nu * nu;
            Complex t = Complex.Log(2d / z), mu = nu * t;
            (ddouble g1, ddouble g2) = Gamma12(nu);

            Complex f = (g1 * Complex.Cosh(mu) + g2 * t * Sinhc(mu)) / ddouble.Sinc(nu);
            Complex r = Complex.Pow(z / 2d, nu);
            Complex p = ddouble.Gamma(1d + nu) / (r * 2d), q = ddouble.Gamma(1d - nu) * r * 0.5d;

            Complex c = f, v = Complex.Ldexp(z * z, -2), u = v;

            for (int k = 1; k <= terms; k++) {
                f = (k * f + p + q) / (k * k - nu2);
                c = SeriesUtil.Add(c, u, f, out bool convergence);

                if (convergence) {
                    break;
                }

                p /= k - nu;
                q /= k + nu;
                u *= v / (k + 1);

                IterationLogger.Log("AmosPowerSeries", k);
            }

            return c;
        }

        public static (ddouble g1, ddouble g2) Gamma12(ddouble nu) {
            Debug.Assert(ddouble.Abs(nu) <= 0.5);

            ddouble nu2 = nu * nu;
            ddouble gm = 1 / ddouble.Gamma(1d - nu), gp = 1 / ddouble.Gamma(1d + nu);

            ddouble g1 = G1Coef[^1];

            for (int i = G1Coef.Count - 2; i >= 0; i--) {
                g1 = g1 * nu2 + G1Coef[i];
            }

            ddouble g2 = (gm + gp) / 2d;

            return (g1, g2);
        }

        public static Complex Sinhc(Complex z) {
            if (ddouble.ILogB(z.Magnitude) > -53) {
                return Complex.Sinh(z) / z;
            }
            else {
                Complex z2 = z * z;

                return (120d + z2 * (20d + z2)) / 120d;
            }
        }

        public static readonly ReadOnlyCollection<ddouble> G1Coef = new([
            (-1, -1, 0x93C467E37DB0C7A4uL, 0xD1BE3F810152CB57uL),
            (+1, -5, 0xAC0AF47D13823E47uL, 0xA15A643C9C84B042uL),
            (+1, -5, 0xACD7881E1A0493DFuL, 0x1B4048388CACA42FuL),
            (-1, -8, 0xEC8CE293FB058CADuL, 0xD1D69DD3DC1E1318uL),
            (+1, -13, 0xE1B27F378AB1E74CuL, 0x7C04703C8691B6EDuL),
            (+1, -16, 0xA8E7457A3F55EFEDuL, 0x92173AFE34B827DBuL),
            (-1, -20, 0x981284EDE06F1640uL, 0xF6FFCB7E4E64F0B8uL),
            (-1, -28, 0xD225BDD116B14565uL, 0x21CC6FD93419DC08uL),
            (+1, -30, 0xA25A676E51C47BE3uL, 0x89819C393897FFDEuL),
            (-1, -37, 0x88E832DFD7833A2DuL, 0x6B19F0B77E71C666uL),
            (-1, -41, 0x8F900A8991E681C8uL, 0xF6862A8BDDBA9233uL),
            (+1, -48, 0xC0AFB9502F8B79A4uL, 0xCB318673047D37D7uL),
            (+1, -53, 0x88309EEF2BD445EAuL, 0x75AE6E0286820762uL),
            (-1, -60, 0xD06E3B87DC524CDAuL, 0x47E432D83C3B6F75uL),
            (-1, -66, 0xA1EBCD25C8674023uL, 0x9930EA5F35E39803uL),
            (+1, -73, 0xF85772F7D97E5EBEuL, 0xEE622174D385FF99uL),
        ]);
    }
}
