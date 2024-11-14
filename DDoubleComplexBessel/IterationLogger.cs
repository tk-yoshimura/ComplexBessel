using System.Collections.Concurrent;
using System.Diagnostics;

namespace DDoubleComplexBessel {
    static class IterationLogger {
        static readonly ConcurrentDictionary<string, int> table = [];

        public static void Log(string tag, int n) {
            if (!table.TryGetValue(tag, out int k)) {
                table[tag] = n;

                Trace.WriteLine($"{tag}: {n}");
            }
            else if (k < n) {
                table[tag] = n;

                Trace.WriteLine($"{tag}: {n}");
            }
        }
    }
}
