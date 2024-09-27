using ComplexBessel;
using MultiPrecision;
using MultiPrecisionComplex;

namespace ComplexBesselTests {
    [TestClass()]
    public class BesselItoJTests {
        
        readonly Complex<Pow2.N4>[] zs = [
            (1, 0), (1, 0.5), (2, 0), (2, 0.5), (2, 1), (0, 1), (0.5, 1), (0, 2), (0.5, 2), (1, 2)
        ];

        [TestMethod()]
        public void BesselItoJTest() {
            for (double nu = -8; nu <= 8; nu += 0.25) {
                Console.WriteLine(nu);

                Complex<Pow2.N4> c = (MultiPrecision<Pow2.N4>.CosPI(nu / 2), -MultiPrecision<Pow2.N4>.SinPI(nu / 2));

                foreach (Complex<Pow2.N4> z in zs) {
                    Complex<Pow2.N4> bj = PowerSeries<Pow2.N4>.BesselJ(nu, z);
                    Complex<Pow2.N4> bi = PowerSeries<Pow2.N4>.BesselI(nu, (z.I, z.R));

                    Complex<Pow2.N4> y = (c * bi).Conj;

                    Console.WriteLine($"{z}\n {bj}\n {y}");

                    Assert.IsTrue((y - bj).Magnitude / bj.Magnitude < 1e-30);
                }

                Console.WriteLine(string.Empty);
            }
        }
    }
}