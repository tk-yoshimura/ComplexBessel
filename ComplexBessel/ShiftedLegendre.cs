using System.Collections.Concurrent;
using System.Numerics;

namespace ComplexBessel {
    class ShiftedLegendre {
        private static readonly ConcurrentDictionary<(int m, int n), BigInteger> table = [];

        static ShiftedLegendre() {
            table[(0, 0)] = 1;
        }

        public static BigInteger Table(int n, int k) {
            if (n < 0 || k < 0 || n < k) {
                throw new ArgumentOutOfRangeException($"{nameof(n)},{nameof(k)}", $"{nameof(n)}>={nameof(k)}>=0");
            }

            if (!table.TryGetValue((n, k), out BigInteger v)) {
                v = -((n == k) ? (Table(n, k - 1) * 2 / n) : (Table(n - 1, k) * (n + k) / (n - k)));
                table[(n, k)] = v;
            }

            return v;
        }
    }
}
