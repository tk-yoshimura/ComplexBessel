﻿using MultiPrecision;

namespace ComplexBesselBenchmark {
    internal struct Plus1<N> : IConstant where N : struct, IConstant {
        public readonly int Value => checked(default(N).Value + 1);
    }
}
