﻿using DoubleDouble;
using DoubleDoubleComplex;
using System.Runtime.CompilerServices;

namespace DDoubleOptimizedBessel {
    internal static class SeriesUtil {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ddouble Add(ddouble c, ddouble s, ddouble a, out bool convergence) {
            ddouble c_next = c + s * a;

            convergence = c == c_next;

            return c_next;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ddouble Add(ddouble c, ddouble s, ddouble a, ddouble b, out bool convergence) {
            ddouble c_next_x = c + s * (a + b), c_next_y = c + s * (a - b);

            convergence = (c == c_next_x) && (c == c_next_y);

            return c_next_x;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Complex Add(Complex c, Complex s, Complex a, out bool convergence) {
            Complex c_next = c + s * a;

            convergence = c == c_next;

            return c_next;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Complex Add(Complex c, Complex s, Complex a, Complex b, out bool convergence) {
            Complex c_next_x = c + s * (a + b), c_next_y = c + s * (a - b);

            convergence = (c == c_next_x) && (c == c_next_y);

            return c_next_x;
        }
    }
}
