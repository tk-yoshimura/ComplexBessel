using MultiPrecision;
using MultiPrecisionComplex;
using System.Runtime.CompilerServices;

namespace ComplexBessel {
    internal static class SeriesUtil<N> where N : struct, IConstant {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Complex<N> Add(Complex<N> c, Complex<N> s, Complex<N> a, out bool convergence) {
            Complex<N> c_next = c + s * a;

            convergence = c == c_next;

            return c_next;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Complex<N> Add(Complex<N> c, Complex<N> s, Complex<N> a, Complex<N> b, out bool convergence) {
            Complex<N> c_next_x = c + s * (a + b), c_next_y = c + s * (a - b);

            convergence = (c == c_next_x) && (c == c_next_y);

            return c_next_x;
        }
    }
}
